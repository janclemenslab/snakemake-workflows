"""Workflow-aware helpers for dashboard-based Snakemake submission."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re

from fabric import Connection

from . import project


_VIDEO_SUFFIXES = frozenset({".mp4", ".avi"})
_RULE_RE = re.compile(r"^(?:localrule|rule|checkpoint)\s+([^:]+):$", re.MULTILINE)
_LOG_PATH_RE = re.compile(r"^Log:\s*(.+)$", re.MULTILINE)
_PID_RE = re.compile(r"^Started Snakemake controller PID\s+(\d+)", re.MULTILINE)


@dataclass(frozen=True)
class SubmissionLaunchResult:
    ok: bool
    message: str
    detail: str = ""
    controller_log: str = ""
    controller_pid: str = ""
    stdout: str = ""
    stderr: str = ""
    target_count: int = 0


def normalize_submission_targets(*chunks: str) -> list[str]:
    targets: list[str] = []
    seen: set[str] = set()
    for chunk in chunks:
        for raw_line in chunk.splitlines():
            target = raw_line.strip()
            if not target or target in seen:
                continue
            seen.add(target)
            targets.append(target)
    return targets


def submission_metadata(project_dir: Path) -> dict[str, object]:
    dat_dir = project_dir / "dat"
    snakefile = _workflow_snakefile(project_dir)
    dat_mtime_ns = dat_dir.stat().st_mtime_ns if dat_dir.exists() else -1
    snakefile_mtime_ns = snakefile.stat().st_mtime_ns if snakefile and snakefile.exists() else -1
    return _submission_metadata_cached(
        str(project_dir.resolve()),
        dat_mtime_ns,
        snakefile_mtime_ns,
    )


@lru_cache(maxsize=8)
def _submission_metadata_cached(
    project_dir_str: str,
    dat_mtime_ns: int,
    snakefile_mtime_ns: int,
) -> dict[str, object]:
    del dat_mtime_ns, snakefile_mtime_ns
    project_dir = Path(project_dir_str)

    if project_dir.name == "playback":
        return _playback_submission_metadata(project_dir)
    if project_dir.name == "chainingmic":
        return _chainingmic_submission_metadata(project_dir)
    return _generic_submission_metadata(project_dir)


def launch_dashboard_submission(
    *,
    project_dir: Path,
    user: str,
    host: str,
    targets: list[str] | None = None,
) -> SubmissionLaunchResult:
    queue_user = (user or "").strip()
    queue_host = (host or "").strip() or _default_remote_host()

    if not queue_user:
        return SubmissionLaunchResult(
            ok=False,
            message="Set the SSH user before submitting.",
        )
    if not queue_host:
        return SubmissionLaunchResult(
            ok=False,
            message="Set the SSH host before submitting.",
        )

    project.configure_project(project_dir.name, project_dir=project_dir)
    controller_cmd = project.snakemake_submit_command(
        targets=targets or [],
        force=bool(targets),
    )
    launch_cmd = project.controller_launch_command(controller_cmd)

    try:
        with Connection(host=queue_host, user=queue_user) as connection:
            result = connection.run(launch_cmd, hide=True, warn=True)
    except Exception as exc:
        return SubmissionLaunchResult(
            ok=False,
            message="Remote submission failed before the controller could start.",
            detail=str(exc),
            target_count=len(targets or []),
        )

    stdout = (result.stdout or "").strip()
    stderr = (result.stderr or "").strip()
    log_path = _normalize_controller_log_path(project_dir, _extract_match(_LOG_PATH_RE, stdout))
    controller_pid = _extract_match(_PID_RE, stdout)

    if not result.ok:
        return SubmissionLaunchResult(
            ok=False,
            message="Snakemake submission failed.",
            detail=stderr or stdout or "The remote command exited with a non-zero status.",
            controller_log=log_path,
            controller_pid=controller_pid,
            stdout=stdout,
            stderr=stderr,
            target_count=len(targets or []),
        )

    if targets:
        message = f"Submitted {len(targets)} explicit target{'s' if len(targets) != 1 else ''}."
    else:
        message = "Submitted the workflow default target."

    detail_parts = []
    if controller_pid:
        detail_parts.append(f"Controller PID {controller_pid}")
    if log_path:
        detail_parts.append(Path(log_path).name)

    return SubmissionLaunchResult(
        ok=True,
        message=message,
        detail=" | ".join(detail_parts),
        controller_log=log_path,
        controller_pid=controller_pid,
        stdout=stdout,
        stderr=stderr,
        target_count=len(targets or []),
    )


def _default_remote_host() -> str:
    if project.DEFAULT_HOSTS:
        return (project.DEFAULT_HOSTS[0].get("host") or "").strip()
    return ""


def _workflow_snakefile(project_dir: Path) -> Path | None:
    for candidate in (
        project_dir / "workflow" / "snakefile",
        project_dir / "Snakefile",
        project_dir / "snakefile",
    ):
        if candidate.exists():
            return candidate
    return None


def _generic_submission_metadata(project_dir: Path) -> dict[str, object]:
    snakefile = _workflow_snakefile(project_dir)
    rules = []
    if snakefile is not None:
        for rule_name in sorted(set(_RULE_RE.findall(snakefile.read_text())), key=str.casefold):
            rules.append(
                {
                    "name": rule_name,
                    "label": rule_name,
                    "description": "Manual target only",
                }
            )

    return {
        "project_name": project_dir.name,
        "kind": "manual",
        "rules": rules,
        "experiments": [],
        "quick_submit_help": "Submit everything that is currently pending in the workflow.",
        "selective_submit_help": (
            "This workflow does not expose curated file mappings yet. "
            "Add curated experiment mappings if you want checklist-based submission."
        ),
    }


def _playback_submission_metadata(project_dir: Path) -> dict[str, object]:
    rules = [
        {
            "name": "detect_fly_chambers",
            "label": "Detect chambers",
            "description": "Create chamber detections for each selected video.",
        },
        {
            "name": "merge_tracks",
            "label": "Track + merge",
            "description": "Build merged SLEAP tracks for each selected video.",
        },
        {
            "name": "pb_speed",
            "label": "Speed traces",
            "description": "Build speed traces and playlists for each selected video.",
        },
        {
            "name": "report_plots",
            "label": "Report pages",
            "description": "Build per-video HTML report pages.",
        },
    ]

    experiments = []
    for session_name, source_files in _video_groups(project_dir):
        for video_name, source_path in source_files:
            targets = {
                "detect_fly_chambers": [
                    f"dat/{session_name}/{video_name}_fly_chambers/{video_name}_detections.json"
                ],
                "merge_tracks": [f"res/{session_name}/{video_name}_sleap.h5"],
                "pb_speed": [f"res/{session_name}/{video_name}_spd.npz"],
                "report_plots": [
                    f"res/{session_name}/{video_name}_report/{video_name}_report.html"
                ],
            }
            experiments.append(
                {
                    "id": f"experiment:{session_name}/{video_name}",
                    "label": video_name,
                    "group": session_name,
                    "path": source_path,
                    "targets": targets,
                    "rule_status": {
                        rule_name: _targets_exist(project_dir, rule_targets)
                        for rule_name, rule_targets in targets.items()
                    },
                }
            )

    return {
        "project_name": project_dir.name,
        "kind": "curated",
        "rules": rules,
        "experiments": experiments,
        "quick_submit_help": "Equivalent to running the workflow default target remotely.",
        "selective_submit_help": (
            "Choose processing stages and one or more sessions/videos. "
            "Finished experiments remain selectable so you can rerun them deliberately."
        ),
    }


def _chainingmic_submission_metadata(project_dir: Path) -> dict[str, object]:
    rules = [
        {
            "name": "sleap",
            "label": "SLEAP",
            "description": "Build tracking output for each selected recording.",
        },
        {
            "name": "song",
            "label": "Song",
            "description": "Build DAS annotations for each selected recording.",
        },
    ]

    experiments = []
    for folder_name, source_files in _video_groups(project_dir):
        for sample_name, source_path in source_files:
            targets = {
                "sleap": [f"res/{folder_name}/{sample_name}_sleap.h5"],
                "song": [f"res/{folder_name}/{sample_name}_annotations.csv"],
            }
            experiments.append(
                {
                    "id": f"experiment:{folder_name}/{sample_name}",
                    "label": sample_name,
                    "group": folder_name,
                    "path": source_path,
                    "targets": targets,
                    "rule_status": {
                        rule_name: _targets_exist(project_dir, rule_targets)
                        for rule_name, rule_targets in targets.items()
                    },
                }
            )

    return {
        "project_name": project_dir.name,
        "kind": "curated",
        "rules": rules,
        "experiments": experiments,
        "quick_submit_help": "Equivalent to running the workflow default target remotely.",
        "selective_submit_help": (
            "Choose stages plus experiments. "
            "Finished experiments remain selectable so you can rerun them deliberately."
        ),
    }


def _video_groups(project_dir: Path) -> list[tuple[str, list[tuple[str, str]]]]:
    dat_dir = project_dir / "dat"
    if not dat_dir.exists():
        return []

    groups = []
    for folder in sorted(
        (path for path in dat_dir.iterdir() if path.is_dir() and not path.name.startswith(".")),
        key=lambda path: path.name.casefold(),
    ):
        source_by_name: dict[str, Path] = {}
        prefix = f"{folder.name}_"
        for candidate in folder.iterdir():
            if not candidate.is_file() or candidate.name.startswith("."):
                continue
            if candidate.suffix.lower() not in _VIDEO_SUFFIXES:
                continue

            sample_name = candidate.stem
            if sample_name != folder.name and not sample_name.startswith(prefix):
                continue

            current = source_by_name.get(sample_name)
            if current is None or (
                current.suffix.lower() != ".mp4" and candidate.suffix.lower() == ".mp4"
            ):
                source_by_name[sample_name] = candidate

        if not source_by_name:
            continue

        groups.append(
            (
                folder.name,
                [
                    (
                        sample_name,
                        _relative_path(source_path, project_dir),
                    )
                    for sample_name, source_path in sorted(
                        source_by_name.items(),
                        key=lambda entry: entry[0].casefold(),
                    )
                ],
            )
        )
    return groups


def _relative_path(path: Path, project_dir: Path) -> str:
    return str(path.relative_to(project_dir)).replace("\\", "/")


def _targets_exist(project_dir: Path, targets: list[str]) -> bool:
    if not targets:
        return False
    return all((project_dir / target).exists() for target in targets)


def _extract_match(pattern: re.Pattern[str], text: str) -> str:
    match = pattern.search(text or "")
    return match.group(1).strip() if match else ""


def _normalize_controller_log_path(project_dir: Path, log_path: str) -> str:
    if not log_path:
        return ""
    candidate = Path(log_path)
    if candidate.is_absolute():
        return str(candidate)
    return str((project_dir / candidate).resolve())
