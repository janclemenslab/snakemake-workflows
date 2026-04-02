import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import av
import numpy as np
from PIL import Image, ImageDraw, ImageFilter


@dataclass
class Component:
    area: int
    bbox: tuple[int, int, int, int]
    centroid: tuple[float, float]
    max_score: float


@dataclass
class ChamberDetection:
    index: int
    bbox: tuple[int, int, int, int]
    detection_bbox: tuple[int, int, int, int]
    has_fly: bool
    fly_component: Component | None


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values.astype(np.float32, copy=False)
    kernel = np.ones(window, dtype=np.float32) / float(window)
    return np.convolve(values.astype(np.float32), kernel, mode="same")


def threshold_from_profile(profile: np.ndarray, fraction: float) -> float:
    minimum = float(np.min(profile))
    maximum = float(np.max(profile))
    return minimum + fraction * (maximum - minimum)


def find_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(mask.tolist()):
        if value and start is None:
            start = idx
        elif not value and start is not None:
            runs.append((start, idx))
            start = None
    if start is not None:
        runs.append((start, int(mask.size)))
    return runs


def largest_run(mask: np.ndarray) -> tuple[int, int]:
    runs = find_runs(mask)
    if not runs:
        raise RuntimeError("No bright plate region found in frame profile.")
    return max(runs, key=lambda run: run[1] - run[0])


def read_first_frame(video_path: Path) -> np.ndarray:
    with av.open(str(video_path)) as container:
        frame = next(container.decode(video=0))
    return frame.to_ndarray(format="rgb24")


def to_grayscale(frame_rgb: np.ndarray) -> np.ndarray:
    frame = frame_rgb.astype(np.float32)
    return 0.299 * frame[:, :, 0] + 0.587 * frame[:, :, 1] + 0.114 * frame[:, :, 2]


def clamp_bbox(
    bbox: tuple[int, int, int, int], width: int, height: int
) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = bbox
    return (
        max(0, min(x0, width - 1)),
        max(0, min(y0, height - 1)),
        max(1, min(x1, width)),
        max(1, min(y1, height)),
    )


def expand_bbox(
    bbox: tuple[int, int, int, int], padding: int, width: int, height: int
) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = bbox
    return clamp_bbox(
        (x0 - padding, y0 - padding, x1 + padding, y1 + padding), width, height
    )


def resize_grayscale(image: np.ndarray, size: tuple[int, int]) -> np.ndarray:
    return np.asarray(
        Image.fromarray(np.clip(image, 0, 255).astype(np.uint8)).resize(
            size, resample=Image.Resampling.BILINEAR
        ),
        dtype=np.float32,
    )


def max_window_sum_bbox(
    values: np.ndarray, window_height: int, window_width: int
) -> tuple[float, tuple[int, int, int, int]]:
    height, width = values.shape
    if height < window_height or width < window_width:
        return float(values.sum()), (0, 0, width, height)

    integral = np.pad(values, ((1, 0), (1, 0)), mode="constant").cumsum(0).cumsum(1)
    sums = (
        integral[window_height:, window_width:]
        - integral[:-window_height, window_width:]
        - integral[window_height:, :-window_width]
        + integral[:-window_height, :-window_width]
    )
    top, left = np.unravel_index(np.argmax(sums), sums.shape)
    return float(sums[top, left]), (
        left,
        top,
        left + window_width,
        top + window_height,
    )


def detect_plate_bbox(gray: np.ndarray) -> tuple[int, int, int, int]:
    row_profile = moving_average(gray.mean(axis=1), 31)
    row_threshold = threshold_from_profile(row_profile, 0.55)
    y0, y1 = largest_run(row_profile > row_threshold)

    cropped = gray[y0:y1, :]
    col_profile = moving_average(cropped.mean(axis=0), 201)
    col_threshold = threshold_from_profile(col_profile, 0.35)
    x0, x1 = largest_run(col_profile > col_threshold)

    pad_x = max(12, int((x1 - x0) * 0.01))
    pad_y = max(12, int((y1 - y0) * 0.02))
    return clamp_bbox(
        (x0 - pad_x, y0 - pad_y, x1 + pad_x, y1 + pad_y),
        gray.shape[1],
        gray.shape[0],
    )


def detect_chambers(
    gray: np.ndarray, plate_bbox: tuple[int, int, int, int]
) -> list[tuple[int, int, int, int]]:
    x0, y0, x1, y1 = plate_bbox
    plate = gray[y0:y1, x0:x1]
    height, width = plate.shape

    core_top = int(height * 0.15)
    core_bottom = int(height * 0.85)
    column_profile = moving_average(plate[core_top:core_bottom, :].mean(axis=0), 9)
    column_threshold = threshold_from_profile(column_profile, 0.58)
    column_runs = find_runs(column_profile > column_threshold)

    min_width = max(20, int(width * 0.03))
    chamber_runs = [run for run in column_runs if run[1] - run[0] >= min_width]
    chamber_specs: list[tuple[int, int, int, int]] = []

    for start, end in chamber_runs:
        chamber_slice = plate[:, start:end]
        row_profile = moving_average(chamber_slice.mean(axis=1), 15)
        row_threshold = threshold_from_profile(row_profile, 0.58)
        try:
            chamber_y0, chamber_y1 = largest_run(row_profile > row_threshold)
        except RuntimeError:
            continue

        chamber_specs.append((start, end, chamber_y0, chamber_y1))

    if not chamber_specs:
        return []

    # Chambers in a plate share a common vertical extent. Using the median
    # top/bottom avoids one-off truncation when a single chamber profile is
    # interrupted by a fly near the edge or uneven illumination.
    common_y0 = int(round(np.median([spec[2] for spec in chamber_specs])))
    common_y1 = int(round(np.median([spec[3] for spec in chamber_specs])))
    pad_y = max(4, int((common_y1 - common_y0) * 0.02))

    chambers: list[tuple[int, int, int, int]] = []
    for start, end, _, _ in chamber_specs:
        pad_x = max(4, int((end - start) * 0.05))
        chambers.append(
            clamp_bbox(
                (
                    x0 + start - pad_x,
                    y0 + common_y0 - pad_y,
                    x0 + end + pad_x,
                    y0 + common_y1 + pad_y,
                ),
                gray.shape[1],
                gray.shape[0],
            )
        )

    return chambers


def detect_fly_component(chamber_gray: np.ndarray) -> Component | None:
    chamber_height, chamber_width = chamber_gray.shape
    scaled_width = max(8, chamber_width // 2)
    scaled_height = max(8, chamber_height // 2)
    scaled = resize_grayscale(chamber_gray, (scaled_width, scaled_height))

    # Downsampling suppresses the chamber mesh, then the local darkness response
    # asks whether there is a fly-sized cluster of dark pixels anywhere inside.
    blurred = np.asarray(
        Image.fromarray(np.clip(scaled, 0, 255).astype(np.uint8)).filter(
            ImageFilter.GaussianBlur(radius=6)
        ),
        dtype=np.float32,
    )
    darkness = np.maximum(0.0, blurred - scaled)

    border = 1
    inner = darkness[border:-border, border:-border]
    threshold = max(3.0, float(np.mean(inner) + 1.25 * np.std(inner)))
    response = np.maximum(0.0, darkness - threshold)
    response[:border, :] = 0.0
    response[-border:, :] = 0.0
    response[:, :border] = 0.0
    response[:, -border:] = 0.0

    score, bbox_small = max_window_sum_bbox(response, window_height=16, window_width=10)
    if score < 1227.0:
        return None

    scale_x = chamber_width / float(scaled_width)
    scale_y = chamber_height / float(scaled_height)
    sx0, sy0, sx1, sy1 = bbox_small
    bbox = (
        max(0, int(math.floor(sx0 * scale_x))),
        max(0, int(math.floor(sy0 * scale_y))),
        min(chamber_width, int(math.ceil(sx1 * scale_x))),
        min(chamber_height, int(math.ceil(sy1 * scale_y))),
    )
    centroid = (
        0.5 * (bbox[0] + bbox[2]),
        0.5 * (bbox[1] + bbox[3]),
    )
    return Component(
        area=max(1, (bbox[2] - bbox[0]) * (bbox[3] - bbox[1])),
        bbox=bbox,
        centroid=centroid,
        max_score=score,
    )


def classify_chambers(
    gray: np.ndarray,
    chambers: Iterable[tuple[int, int, int, int]],
    detection_padding: int,
) -> list[ChamberDetection]:
    detections: list[ChamberDetection] = []
    frame_height, frame_width = gray.shape
    for index, bbox in enumerate(chambers, start=1):
        detection_bbox = expand_bbox(bbox, detection_padding, frame_width, frame_height)
        x0, y0, x1, y1 = detection_bbox
        chamber = gray[y0:y1, x0:x1]
        component = detect_fly_component(chamber)
        detections.append(
            ChamberDetection(
                index=index,
                bbox=bbox,
                detection_bbox=detection_bbox,
                has_fly=component is not None,
                fly_component=component,
            )
        )
    return detections


def draw_debug(
    frame_rgb: np.ndarray,
    plate_bbox: tuple[int, int, int, int],
    detections: Iterable[ChamberDetection],
    output_path: Path,
) -> None:
    image = Image.fromarray(frame_rgb)
    draw = ImageDraw.Draw(image)
    draw.rectangle(plate_bbox, outline=(255, 215, 0), width=4)

    for detection in detections:
        chamber_color = (80, 220, 120) if detection.has_fly else (220, 80, 80)
        draw.rectangle(detection.bbox, outline=chamber_color, width=3)
        label_anchor = (detection.bbox[0] + 4, max(0, detection.bbox[1] - 18))
        draw.text(label_anchor, f"{detection.index}", fill=chamber_color)

        if detection.fly_component is not None:
            fx0, fy0, fx1, fy1 = detection.fly_component.bbox
            cx0, cy0, _, _ = detection.detection_bbox
            draw.rectangle(
                (cx0 + fx0, cy0 + fy0, cx0 + fx1, cy0 + fy1),
                outline=(30, 144, 255),
                width=2,
            )

    image.save(output_path)


def crop_video_region(
    video_path: Path,
    output_path: Path,
    bbox: tuple[int, int, int, int],
    force: bool,
) -> None:
    if output_path.exists() and not force:
        return

    x0, y0, x1, y1 = bbox
    width = x1 - x0
    height = y1 - y0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "ffmpeg",
        "-y" if force else "-n",
        "-i",
        str(video_path),
        "-vf",
        f"crop={width}:{height}:{x0}:{y0}",
        "-an",
        "-c:v",
        "libx264",
        "-crf",
        "18",
        "-preset",
        "fast",
        str(output_path),
    ]
    subprocess.run(
        cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )


def chamber_record_from_manifest(manifest: dict, chamber_index: int) -> dict:
    for chamber in manifest.get("chambers", []):
        if int(chamber.get("index", -1)) == chamber_index:
            return chamber
    raise ValueError(f"Could not find chamber {chamber_index:02d} in detections manifest.")


def get_param(name, default=None):
    try:
        return getattr(snakemake.params, name)  # noqa: F821
    except Exception:
        return default


def main():
    video_path = Path(str(snakemake.input.video)).resolve()  # noqa: F821
    mode = str(get_param("mode", "detect"))
    min_chambers = int(get_param("min_chambers", 6))
    crop_padding = int(get_param("crop_padding", 10))
    detection_padding = int(get_param("detection_padding", 18))
    force = bool(get_param("force", False))
    skip_cropping = bool(get_param("skip_cropping", False))

    if mode == "crop":
        manifest_path = Path(str(snakemake.input.manifest)).resolve()  # noqa: F821
        output_video = Path(str(snakemake.output[0])).resolve()  # noqa: F821
        chamber_index = int(str(snakemake.wildcards.chamber))  # noqa: F821

        with manifest_path.open() as handle:
            manifest = json.load(handle)

        chamber = chamber_record_from_manifest(manifest, chamber_index)
        if not chamber.get("has_fly"):
            raise ValueError(f"Chamber {chamber_index:02d} is not fly-positive in {manifest_path}.")

        crop_bbox = tuple(int(value) for value in chamber["crop_bbox_xyxy"])
        crop_video_region(video_path, output_video, crop_bbox, force=force)
        return

    if mode != "detect":
        raise ValueError(f"Unsupported split_videos mode: {mode}")

    manifest_path = Path(str(snakemake.output.manifest)).resolve()  # noqa: F821
    debug_path = Path(str(snakemake.output.detections_png)).resolve()  # noqa: F821
    output_dir = manifest_path.parent

    if output_dir.parent != video_path.parent:
        raise ValueError(
            f"Expected output dir next to input video. Got {output_dir} for {video_path}."
        )

    frame_rgb = read_first_frame(video_path)
    frame_height, frame_width = frame_rgb.shape[:2]
    gray = to_grayscale(frame_rgb)
    plate_bbox = detect_plate_bbox(gray)
    chambers = detect_chambers(gray, plate_bbox)
    if len(chambers) < min_chambers:
        raise RuntimeError(
            f"Only detected {len(chambers)} chambers in {video_path}, expected at least "
            f"{min_chambers}."
        )

    detections = classify_chambers(gray, chambers, detection_padding)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "video": str(video_path),
        "frame_shape": list(frame_rgb.shape),
        "plate_bbox_xyxy": list(plate_bbox),
        "num_detected_chambers": len(chambers),
        "num_fly_chambers": int(sum(d.has_fly for d in detections)),
        "chambers": [],
    }

    for detection in detections:
        crop_bbox = expand_bbox(detection.bbox, crop_padding, frame_width, frame_height)
        chamber_record = {
            "index": detection.index,
            "bbox_xyxy": list(detection.bbox),
            "detection_bbox_xyxy": list(detection.detection_bbox),
            "crop_bbox_xyxy": list(crop_bbox),
            "has_fly": detection.has_fly,
        }
        if detection.fly_component is not None:
            chamber_record["fly_component"] = {
                "bbox_xyxy": list(detection.fly_component.bbox),
                "area": detection.fly_component.area,
                "centroid_xy": list(detection.fly_component.centroid),
                "max_score": detection.fly_component.max_score,
            }
        manifest["chambers"].append(chamber_record)

        if detection.has_fly and not skip_cropping:
            output_video = (
                output_dir / f"{video_path.stem}_chamber{detection.index:02d}.mp4"
            )
            crop_video_region(video_path, output_video, crop_bbox, force=force)

    output_dir.mkdir(parents=True, exist_ok=True)
    draw_debug(frame_rgb, plate_bbox, detections, debug_path)
    manifest_path.write_text(json.dumps(manifest, indent=2))


main()
