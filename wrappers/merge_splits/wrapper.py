import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np


CHAMBER_PATTERN = re.compile(r"_chamber(\d+)_sleap\.h5$")


def chamber_index_from_path(path: Path) -> int:
    match = CHAMBER_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"Could not infer chamber index from {path}")
    return int(match.group(1))


def decoded_scalar(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return value.decode("utf-8")
    return value


def track_names_for_chamber(
    chamber_index: int,
    raw_track_names: np.ndarray | None,
    track_count: int,
) -> list[str]:
    names = []
    for local_index in range(track_count):
        original_name = ""
        if raw_track_names is not None and local_index < raw_track_names.shape[0]:
            original_name = str(decoded_scalar(raw_track_names[local_index])).strip()
        if original_name:
            names.append(f"chamber{chamber_index:02d}_{original_name}")
        else:
            names.append(f"chamber{chamber_index:02d}_track{local_index:02d}")
    return names


def bytes_array(values: list[str]) -> np.ndarray:
    if not values:
        return np.asarray([], dtype="S1")
    width = max(1, max(len(value.encode("utf-8")) for value in values))
    return np.asarray([value.encode("utf-8") for value in values], dtype=f"S{width}")


def validate_reference_dataset(
    reference_name: str,
    reference_value: np.ndarray,
    candidate_value: np.ndarray,
    path: Path,
) -> None:
    if reference_value.shape != candidate_value.shape or not np.array_equal(
        reference_value, candidate_value
    ):
        raise ValueError(
            f"Dataset {reference_name} in {path} does not match the first chamber track file."
        )


def write_string_dataset(handle: h5py.File, name: str, value: str) -> None:
    handle.create_dataset(name, data=value, dtype=h5py.string_dtype(encoding="utf-8"))


def main():
    output_path = Path(str(snakemake.output.merged)).resolve()  # noqa: F821
    manifest_path = Path(str(snakemake.input.manifest)).resolve()  # noqa: F821
    source_video_path = str(snakemake.input.video)  # noqa: F821
    input_paths = [Path(path).resolve() for path in snakemake.input.chamber_tracks]  # noqa: F821

    if not input_paths:
        raise ValueError("merge_splits requires at least one chamber track input.")

    with manifest_path.open() as handle:
        manifest = json.load(handle)
    chambers_by_index = {int(chamber["index"]): chamber for chamber in manifest["chambers"]}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    reference_arrays: dict[str, np.ndarray] = {}
    merged_tracks = []
    merged_point_scores = []
    merged_instance_scores = []
    merged_tracking_scores = []
    merged_track_occupancy = []
    merged_track_names: list[str] = []
    track_to_chamber: list[int] = []
    input_provenance: list[dict | str] = []
    frame_count = None

    for path in input_paths:
        if not path.exists():
            raise FileNotFoundError(str(path))

        chamber_index = chamber_index_from_path(path)
        chamber = chambers_by_index.get(chamber_index)
        if chamber is None:
            raise ValueError(f"Missing chamber {chamber_index:02d} in {manifest_path}")
        if not chamber.get("has_fly"):
            raise ValueError(f"Chamber {chamber_index:02d} is not marked fly-positive in {manifest_path}")

        crop_x0, crop_y0, _, _ = chamber["crop_bbox_xyxy"]

        with h5py.File(path, "r") as handle:
            for dataset_name in ("edge_inds", "edge_names", "node_names"):
                dataset_value = np.asarray(handle[dataset_name][...])
                if dataset_name not in reference_arrays:
                    reference_arrays[dataset_name] = dataset_value
                else:
                    validate_reference_dataset(
                        dataset_name, reference_arrays[dataset_name], dataset_value, path
                    )

            tracks = np.asarray(handle["tracks"][...])
            point_scores = np.asarray(handle["point_scores"][...])
            instance_scores = np.asarray(handle["instance_scores"][...])
            tracking_scores = np.asarray(handle["tracking_scores"][...])
            track_occupancy = np.asarray(handle["track_occupancy"][...])
            raw_track_names = (
                np.asarray(handle["track_names"][...]) if "track_names" in handle else None
            )

            if tracks.ndim != 4 or tracks.shape[1] != 2:
                raise ValueError(f"Unexpected tracks layout in {path}: {tracks.shape}")
            if point_scores.ndim != 3:
                raise ValueError(f"Unexpected point_scores layout in {path}: {point_scores.shape}")
            if instance_scores.ndim != 2:
                raise ValueError(
                    f"Unexpected instance_scores layout in {path}: {instance_scores.shape}"
                )
            if tracking_scores.ndim != 2:
                raise ValueError(
                    f"Unexpected tracking_scores layout in {path}: {tracking_scores.shape}"
                )
            if track_occupancy.ndim != 2:
                raise ValueError(
                    f"Unexpected track_occupancy layout in {path}: {track_occupancy.shape}"
                )

            local_frame_count = tracks.shape[3]
            if frame_count is None:
                frame_count = local_frame_count
            elif frame_count != local_frame_count:
                raise ValueError(
                    f"Mismatched frame count in {path}: {local_frame_count}, expected {frame_count}"
                )

            track_count = tracks.shape[0]
            if track_occupancy.shape != (frame_count, track_count):
                raise ValueError(
                    f"Unexpected track_occupancy shape in {path}: {track_occupancy.shape}"
                )
            if point_scores.shape[0] != track_count or point_scores.shape[2] != frame_count:
                raise ValueError(
                    f"Unexpected point_scores shape in {path}: {point_scores.shape}"
                )
            if instance_scores.shape != (track_count, frame_count):
                raise ValueError(
                    f"Unexpected instance_scores shape in {path}: {instance_scores.shape}"
                )
            if tracking_scores.shape != (track_count, frame_count):
                raise ValueError(
                    f"Unexpected tracking_scores shape in {path}: {tracking_scores.shape}"
                )

            tracks = tracks.copy()
            tracks[:, 0, :, :] += float(crop_x0)
            tracks[:, 1, :, :] += float(crop_y0)

            merged_tracks.append(tracks)
            merged_point_scores.append(point_scores)
            merged_instance_scores.append(instance_scores)
            merged_tracking_scores.append(tracking_scores)
            merged_track_occupancy.append(track_occupancy)
            merged_track_names.extend(
                track_names_for_chamber(chamber_index, raw_track_names, track_count)
            )
            track_to_chamber.extend([chamber_index] * track_count)

            provenance_value = decoded_scalar(handle["provenance"][()])
            try:
                input_provenance.append(json.loads(provenance_value))
            except Exception:
                input_provenance.append(provenance_value)

    assert frame_count is not None

    merged_tracks_array = np.concatenate(merged_tracks, axis=0)
    merged_point_scores_array = np.concatenate(merged_point_scores, axis=0)
    merged_instance_scores_array = np.concatenate(merged_instance_scores, axis=0)
    merged_tracking_scores_array = np.concatenate(merged_tracking_scores, axis=0)
    merged_track_occupancy_array = np.concatenate(merged_track_occupancy, axis=1)
    merged_track_names_array = bytes_array(merged_track_names)
    track_to_chamber_array = np.asarray(track_to_chamber, dtype=np.int32)

    provenance = {
        "merge_type": "full_frame_materialized",
        "coordinate_space": "full_frame",
        "merge_timestamp": datetime.now(timezone.utc).isoformat(),
        "source_file": source_video_path,
        "manifest_path": str(manifest_path),
        "inputs": [str(path) for path in input_paths],
        "input_provenance": input_provenance,
    }

    with h5py.File(output_path, "w") as out:
        out.attrs["merge_type"] = "full_frame_materialized"
        out.attrs["coordinate_space"] = "full_frame"
        out.attrs["inputs"] = [str(path) for path in input_paths]
        out.attrs["manifest_path"] = str(manifest_path)

        out.create_dataset("edge_inds", data=reference_arrays["edge_inds"])
        out.create_dataset("edge_names", data=reference_arrays["edge_names"])
        out.create_dataset("node_names", data=reference_arrays["node_names"])
        out.create_dataset("instance_scores", data=merged_instance_scores_array)
        write_string_dataset(out, "labels_path", "")
        out.create_dataset("point_scores", data=merged_point_scores_array)
        write_string_dataset(out, "provenance", json.dumps(provenance))
        out.create_dataset("track_names", data=merged_track_names_array)
        out.create_dataset("track_occupancy", data=merged_track_occupancy_array)
        out.create_dataset("tracking_scores", data=merged_tracking_scores_array)
        out.create_dataset("tracks", data=merged_tracks_array)
        out.create_dataset("track_to_chamber", data=track_to_chamber_array)
        out.create_dataset("video_ind", data=np.int64(0))
        write_string_dataset(out, "video_path", source_video_path)


main()
