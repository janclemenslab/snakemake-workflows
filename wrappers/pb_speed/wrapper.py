"""Compute stimulus-aligned speed traces from merged SLEAP tracks."""

import ast
import datetime as dt
import json
import os
import re
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import numpy as np
import pandas as pd
import scipy.signal


TRACK_SUFFIX = "_sleap.h5"
DEFAULTS = {
    "target_sampling_rate": 1_000,
    "resample_video_data": True,
    "speed_median_kernel": 5,
    "song_env_window_seconds": 0.05,
    "song_env_threshold": 1.0,
    "stim_closing_samples": 5_001,
    "prefix_seconds": 5.0,
    "duration_seconds": 10.0,
    "chamber_size_mm": 45.7,
}


def parse_log_file(log_file: Path) -> pd.DataFrame:
    with log_file.open() as handle:
        lines = [line for line in handle.readlines() if "cnt: " in line]

    logs = []
    for line in lines:
        token = line.split(" ")
        log = {
            "timestamp": dt.datetime.strptime(token[0], "%Y-%m-%d,%H:%M:%S.%f"),
            "rig": token[1].rstrip(":"),
        }

        for key, value in zip(token[2:-1:2], token[3:-1:2]):
            key = key.rstrip(":")
            value = value.rstrip(";")
            try:
                parsed_value = ast.literal_eval(value)
            except Exception:
                parsed_value = value
            log[key] = parsed_value

        logs.append(log)

    return pd.DataFrame(logs)


def infer_root(track_path: Path) -> Path:
    for parent in track_path.parents:
        if parent.name == "res":
            return parent.parent
    return track_path.parents[2]


def infer_context(track_path: Path) -> tuple[Path, str, str, Path]:
    if not track_path.name.endswith(TRACK_SUFFIX):
        raise ValueError(
            f"Expected a merged SLEAP track file ending in {TRACK_SUFFIX}: {track_path}"
        )

    root = infer_root(track_path)
    session = track_path.parent.name
    video = track_path.name[: -len(TRACK_SUFFIX)]
    log_path = root / "dat" / session / f"{session}_daq.log"
    return root, session, video, log_path


def chamber_size_px_from_record(chamber: dict) -> float:
    bbox = chamber.get("bbox_xyxy")
    if not isinstance(bbox, (list, tuple)) or len(bbox) != 4:
        raise ValueError(f"Invalid chamber bbox_xyxy entry: {bbox}")

    _, y0, _, y1 = (float(value) for value in bbox)
    chamber_size_px = y1 - y0
    if chamber_size_px <= 0:
        raise ValueError(f"Non-positive chamber height from bbox_xyxy: {bbox}")

    return chamber_size_px


def mean_chamber_size_px(manifest_path: Path) -> float:
    with manifest_path.open() as handle:
        manifest = json.load(handle)

    chambers = manifest.get("chambers")
    if not isinstance(chambers, list) or not chambers:
        raise ValueError(f"No chambers found in detections manifest: {manifest_path}")

    chamber_sizes_px = [chamber_size_px_from_record(chamber) for chamber in chambers]
    return float(np.mean(chamber_sizes_px))


def odd_kernel_size(value: int) -> int:
    value = max(1, int(value))
    if value % 2 == 0:
        value += 1
    return value


def close_small_gaps_1d(x: np.ndarray, max_gap: int) -> np.ndarray:
    """Fill False gaps shorter than max_gap between True runs in a 1D array."""
    x = np.asarray(x, dtype=bool)
    max_gap = int(max_gap)

    if x.size == 0 or max_gap <= 0:
        return x.copy()

    padded = np.concatenate(([False], x, [False]))
    transitions = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(transitions == 1)
    ends = np.flatnonzero(transitions == -1)

    if starts.size <= 1:
        return x.copy()

    gaps = starts[1:] - ends[:-1]
    out = x.copy()
    for gap_idx in np.flatnonzero(gaps < max_gap):
        out[ends[gap_idx] : starts[gap_idx + 1]] = True

    return out


def compute_speed(ds, median_kernel: int) -> np.ndarray:
    import xarray_behave as xb

    positions = np.nanmean(np.asarray(ds.pose_positions_allo), axis=-2)
    velocity = np.asarray(xb.metrics.velocity(positions))
    speed_raw = np.linalg.norm(velocity, axis=-1)

    median_kernel = odd_kernel_size(median_kernel)
    if speed_raw.ndim == 1:
        return scipy.signal.medfilt(speed_raw, kernel_size=median_kernel)
    return scipy.signal.medfilt(speed_raw, kernel_size=(median_kernel, 1))


def detect_stimulus_events(
    ds,
    env_window_seconds: float,
    env_threshold: float,
    stim_closing_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    fs = float(ds.attrs["target_sampling_rate_Hz"])
    win_len = max(1, int(round(env_window_seconds * fs)))
    std = max(1, win_len // 8)

    window = scipy.signal.windows.gaussian(win_len, std)
    window /= np.sum(window)

    song = np.asarray(ds.song_raw[:, 0]).astype(float)
    env = np.sqrt(np.convolve(song**2, window, mode="full"))
    env = env[win_len // 2 : win_len // 2 + song.shape[0]]

    stim = close_small_gaps_1d(env > env_threshold, stim_closing_samples)

    stim_diff = np.diff(stim.astype(float), prepend=0)
    stim_onsets = np.flatnonzero(stim_diff == 1)
    stim_offsets = np.flatnonzero(stim_diff == -1)

    if stim.size and stim[-1]:
        stim_offsets = np.append(stim_offsets, stim.size - 1)

    event_count = min(len(stim_onsets), len(stim_offsets))
    stim_onsets = stim_onsets[:event_count]
    stim_offsets = stim_offsets[:event_count]

    onset_times = np.asarray(ds.sampletime[stim_onsets], dtype=float)
    offset_times = np.asarray(ds.sampletime[stim_offsets], dtype=float)
    return onset_times, offset_times


def extract_aligned_traces(
    speed: np.ndarray,
    time_axis: np.ndarray,
    onset_times: np.ndarray,
    prefix: int,
    duration: int,
) -> np.ndarray:
    speed = np.asarray(speed, dtype=float)
    if speed.ndim == 1:
        speed = speed[:, np.newaxis]

    trace_len = prefix + duration
    traces = np.full((len(onset_times), trace_len, speed.shape[1]), np.nan, dtype=float)

    for idx, onset_time in enumerate(onset_times):
        onset_index = int(np.searchsorted(time_axis, float(onset_time), side="left"))
        start = onset_index - prefix
        stop = onset_index + duration

        src_start = max(start, 0)
        src_stop = min(stop, speed.shape[0])
        if src_stop <= src_start:
            continue

        dst_start = max(0, -start)
        dst_stop = dst_start + (src_stop - src_start)
        traces[idx, dst_start:dst_stop, :] = speed[src_start:src_stop, :]

    return traces


def normalize_log_value(value):
    if isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    if isinstance(value, (pd.Timestamp, dt.datetime, dt.date)):
        return pd.Timestamp(value).isoformat()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value)
    return value


def normalize_dataframe_for_csv(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(lambda column: column.map(normalize_log_value))


def flatten_traces_by_fly(traces: np.ndarray) -> np.ndarray:
    traces = np.asarray(traces, dtype=float)
    return traces.transpose((2, 0, 1)).reshape(-1, traces.shape[1])


def expand_playlist_for_flies(playlist: pd.DataFrame, fly_count: int) -> pd.DataFrame:
    if fly_count <= 0:
        return playlist.iloc[0:0].copy()
    return pd.concat([playlist] * fly_count, ignore_index=True)


def chamber_ids_from_tracks(track_path: Path, fly_count: int) -> list[str]:
    import h5py

    with h5py.File(track_path, "r") as handle:
        raw_names = handle["track_names"][:]

    chamber_ids = []
    for raw_name in raw_names[:fly_count]:
        track_name = (
            raw_name.decode("utf-8") if isinstance(raw_name, bytes) else str(raw_name)
        )
        match = re.search(r"chamber\d+", track_name)
        chamber_ids.append(match.group(0) if match else track_name)

    if len(chamber_ids) != fly_count:
        raise ValueError(
            f"Expected {fly_count} track names in {track_path}, found {len(chamber_ids)}"
        )

    return chamber_ids


# import snakemake
# session = "localhost-20260408_094901"
# video = "localhost-20260408_094901_2"
# track_path = Path(f"res/{session}/{video}_sleap.h5")
# output_path = Path(f"res/{session}/{video}_spd.npz")


def main():
    track_path = Path(str(snakemake.input.tracks)).resolve()  # noqa: F821
    manifest_path = Path(str(snakemake.input.manifest)).resolve()  # noqa: F821
    output_path = Path(str(snakemake.output.speed)).resolve()  # noqa: F821
    playlist_path = Path(str(snakemake.output.playlist)).resolve()  # noqa: F821

    params = DEFAULTS.copy()
    params.update(dict(snakemake.params))  # noqa: F821

    root, session, video, log_path = infer_context(track_path)
    chamber_size_mm = float(params["chamber_size_mm"])
    chamber_size_px = mean_chamber_size_px(manifest_path)
    pixel_size_mm = chamber_size_mm / chamber_size_px

    print(f"Processing merged tracks: {track_path}")
    print(f"Resolved root={root} session={session} video={video}")
    print(f"Using detections manifest: {manifest_path}")
    print(
        f"Using chamber_size_mm={chamber_size_mm:.6f}, "
        f"mean chamber_size_px={chamber_size_px:.6f}, "
        f"pixel_size_mm={pixel_size_mm:.9f}"
    )

    import xarray_behave as xb

    ds = xb.assemble(
        datename=session,
        root=str(root),
        target_sampling_rate=int(params["target_sampling_rate"]),
        resample_video_data=bool(params["resample_video_data"]),
        filepath_poses=str(track_path),
        filepath_video=str(track_path),
        pixel_size_mm=pixel_size_mm,
    )

    speed = compute_speed(ds, median_kernel=int(params["speed_median_kernel"]))
    onset_times, offset_times = detect_stimulus_events(
        ds,
        env_window_seconds=float(params["song_env_window_seconds"]),
        env_threshold=float(params["song_env_threshold"]),
        stim_closing_samples=int(params["stim_closing_samples"]),
    )

    fs = float(ds.attrs["target_sampling_rate_Hz"])
    prefix = int(round(float(params["prefix_seconds"]) * fs))
    duration = int(round(float(params["duration_seconds"]) * fs))
    time_axis = np.asarray(ds.time, dtype=float)
    traces = extract_aligned_traces(
        speed, time_axis, onset_times, prefix=prefix, duration=duration
    )
    fly_count = traces.shape[2]
    chamber_ids = chamber_ids_from_tracks(track_path, fly_count)
    traces = flatten_traces_by_fly(traces)
    rel_time = np.arange(-prefix, duration, dtype=float) / fs

    print(f"Detected {len(onset_times)} stimulus events.")
    print(f"Saving traces with shape {traces.shape} to {output_path}")

    playlist = pd.DataFrame(index=np.arange(len(onset_times)))
    if log_path.exists():
        logs = parse_log_file(log_path).reset_index(drop=True)
        playlist = playlist.join(logs.reindex(playlist.index))
    else:
        print(
            f"DAQ log not found, saving speed traces without log metadata: {log_path}"
        )

    playlist["onset_times"] = onset_times
    playlist["offset_times"] = offset_times
    trial_count = len(playlist)
    playlist = expand_playlist_for_flies(playlist, fly_count)
    playlist.insert(0, "sessionname", session)
    playlist.insert(1, "videoname", video)
    playlist.insert(2, "chamber_id", np.repeat(chamber_ids, trial_count))
    playlist = normalize_dataframe_for_csv(playlist)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    playlist_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Saving playlist with {len(playlist)} rows to {playlist_path}")
    playlist.to_csv(playlist_path, index=False)
    np.savez(
        output_path,
        traces=traces,
        time=rel_time,
        source_tracks=str(track_path),
        session=session,
        video=video,
        chamber_size_mm=chamber_size_mm,
        chamber_size_px=chamber_size_px,
        pixel_size_mm=pixel_size_mm,
    )


try:
    logfile = snakemake.log[0] if snakemake.log else None  # noqa: F821
    if logfile:
        log_path = Path(logfile)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with (
            log_path.open("a") as handle,
            redirect_stdout(handle),
            redirect_stderr(handle),
        ):
            main()
    else:
        main()

except:
    pass
