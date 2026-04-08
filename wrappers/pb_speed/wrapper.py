"""Compute stimulus-aligned speed traces from merged SLEAP tracks."""

from __future__ import annotations

import ast
import datetime as dt
import json
import os
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import numpy as np
import pandas as pd
import scipy.ndimage
import scipy.signal
import xarray_behave as xb


TRACK_SUFFIX = "_sleap.h5"
DEFAULTS = {
    "target_sampling_rate": 1_000,
    "resample_video_data": True,
    "speed_median_kernel": 5,
    "song_env_window_seconds": 0.05,
    "song_env_threshold": 1.0,
    "stim_closing_samples": 1_001,
    "prefix_seconds": 5.0,
    "duration_seconds": 10.0,
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
        raise ValueError(f"Expected a merged SLEAP track file ending in {TRACK_SUFFIX}: {track_path}")

    root = infer_root(track_path)
    session = track_path.parent.name
    video = track_path.name[: -len(TRACK_SUFFIX)]
    log_path = root / "dat" / session / f"{session}_daq.log"
    return root, session, video, log_path


def odd_kernel_size(value: int) -> int:
    value = max(1, int(value))
    if value % 2 == 0:
        value += 1
    return value


def assemble_dataset(
    root: Path,
    session: str,
    track_path: Path,
    target_sampling_rate: int,
    resample_video_data: bool,
):
    return xb.assemble(
        datename=session,
        root=str(root),
        target_sampling_rate=target_sampling_rate,
        resample_video_data=resample_video_data,
        filepath_poses=str(track_path),
        filepath_video=str(track_path),
    )


def compute_speed(ds, median_kernel: int) -> np.ndarray:
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

    stim = scipy.ndimage.binary_closing(
        env > env_threshold,
        np.ones(max(1, int(stim_closing_samples)), dtype=bool),
    )
    stim = np.asarray(stim, dtype=bool)

    stim_diff = np.diff(stim.astype(np.int8), prepend=0)
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


def dataframe_to_npz_kwargs(frame: pd.DataFrame) -> dict[str, np.ndarray]:
    kwargs = {}
    for column in frame.columns:
        values = [normalize_log_value(value) for value in frame[column]]
        kwargs[column] = np.asarray(values)
    return kwargs


def main():
    track_path = Path(str(snakemake.input.tracks)).resolve()  # noqa: F821
    output_path = Path(str(snakemake.output.speed)).resolve()  # noqa: F821

    params = DEFAULTS.copy()
    params.update(dict(snakemake.params))  # noqa: F821

    root, session, video, log_path = infer_context(track_path)
    print(f"Processing merged tracks: {track_path}")
    print(f"Resolved root={root} session={session} video={video}")

    ds = assemble_dataset(
        root=root,
        session=session,
        track_path=track_path,
        target_sampling_rate=int(params["target_sampling_rate"]),
        resample_video_data=bool(params["resample_video_data"]),
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
    traces = extract_aligned_traces(speed, time_axis, onset_times, prefix=prefix, duration=duration)
    rel_time = np.arange(-prefix, duration, dtype=float) / fs

    print(f"Detected {len(onset_times)} stimulus events.")
    print(f"Saving traces with shape {traces.shape} to {output_path}")

    metadata = pd.DataFrame(index=np.arange(len(onset_times)))
    if log_path.exists():
        logs = parse_log_file(log_path).reset_index(drop=True)
        metadata = metadata.join(logs.reindex(metadata.index))
    else:
        print(f"DAQ log not found, saving speed traces without log metadata: {log_path}")

    metadata["onset_times"] = onset_times
    metadata["offset_times"] = offset_times

    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        traces=traces,
        time=rel_time,
        source_tracks=str(track_path),
        session=session,
        video=video,
        **dataframe_to_npz_kwargs(metadata),
    )


logfile = snakemake.log[0] if snakemake.log else None  # noqa: F821
if logfile:
    log_path = Path(logfile)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle, redirect_stdout(handle), redirect_stderr(handle):
        main()
else:
    main()
