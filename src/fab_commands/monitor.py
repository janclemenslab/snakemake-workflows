"""Helpers for monitoring Snakemake controller logs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
import subprocess
from typing import Iterable

from fabric import Connection
from rich.console import Console
from rich.table import Table


_RULE_RE = re.compile(r"^(localrule|rule) ([^:]+):$")
_FIELD_RE = re.compile(r"^\s+([A-Za-z_][\w-]*):\s*(.*)$")
_SUBMIT_RE = re.compile(
    r"Job (\d+) has been submitted with SLURM jobid (\d+) \(log: (.+?)\)\."
)
_FINISH_RE = re.compile(r"Finished jobid: (\d+) \(Rule: ([^)]+)\)")
_PROGRESS_RE = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")
_ERROR_RULE_RE = re.compile(r"^Error in rule ([^:]+):$")
_SLURM_RUN_ID_RE = re.compile(r"^SLURM run ID: (.+)$")
_DISPLAY_LOG_PREFIXES = (
    "INFO:snakemake.logging:",
    "WARNING:snakemake.logging:",
    "ERROR:snakemake.logging:",
)


@dataclass
class JobRecord:
    snakemake_id: str
    rule: str = ""
    kind: str = "rule"
    wildcards: str = ""
    input: str = ""
    output: str = ""
    log: str = ""
    benchmark: str = ""
    slurm_jobid: str = ""
    slurm_log: str = ""
    status: str = "pending"
    slurm_state: str = ""
    elapsed: str = ""
    node: str = ""
    reason: str = ""


@dataclass
class ControllerRun:
    controller_log: Path
    state: str = "unknown"
    slurm_run_id: str = ""
    progress_done: int = 0
    progress_total: int = 0
    progress_percent: int = 0
    status_line: str = ""
    jobs: dict[str, JobRecord] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)

    @property
    def started_at(self) -> str:
        name = self.controller_log.stem
        if not name.startswith("controller-"):
            return self.controller_log.name
        stamp = name.removeprefix("controller-")
        if len(stamp) != 15 or "T" not in stamp:
            return stamp
        return (
            f"{stamp[0:4]}-{stamp[4:6]}-{stamp[6:8]} "
            f"{stamp[9:11]}:{stamp[11:13]}:{stamp[13:15]}"
        )

    def counts(self) -> dict[str, int]:
        counts = {
            "running": 0,
            "failed": 0,
            "finished": 0,
            "submitted": 0,
            "local": 0,
        }
        for job in self.jobs.values():
            if job.kind == "localrule":
                counts["local"] += 1
            if job.status == "failed":
                counts["failed"] += 1
            elif job.status == "finished":
                counts["finished"] += 1
            elif job.status in {"running", "pending", "submitted"}:
                if job.slurm_jobid:
                    counts["submitted"] += 1
                if job.status == "running":
                    counts["running"] += 1
        return counts


def latest_controller_log(project_dir: Path) -> Path | None:
    logs = controller_logs(project_dir)
    if not logs:
        return None
    return logs[-1]


def controller_logs(project_dir: Path) -> list[Path]:
    log_dir = project_dir / "log" / "slurm"
    if not log_dir.exists():
        return []

    logs = sorted(
        path
        for path in log_dir.glob("controller-*.log")
        if not path.name.startswith("._")
    )
    return logs


def _merge_log_values(*values: str) -> str:
    parts: list[str] = []
    seen: set[str] = set()
    for value in values:
        if not value or value == "-":
            continue
        for raw_part in value.split(","):
            cleaned = raw_part.strip()
            cleaned = re.sub(
                r"\s*\(check log file\(s\) for error details\)\s*$",
                "",
                cleaned,
            )
            if not cleaned or cleaned in seen:
                continue
            seen.add(cleaned)
            parts.append(cleaned)
    return ", ".join(parts)


def parse_controller_log(path: Path) -> ControllerRun:
    run = ControllerRun(controller_log=path)
    current_job: dict[str, str] | None = None
    error_mode = False
    error_fields: dict[str, str] = {}
    event_state = "running"
    last_event = "Controller log found."

    def ensure_job(job_id: str) -> JobRecord:
        if job_id not in run.jobs:
            run.jobs[job_id] = JobRecord(snakemake_id=job_id)
        return run.jobs[job_id]

    def finish_error_block() -> None:
        nonlocal error_mode, error_fields, event_state, last_event
        if not error_mode:
            return
        job_id = error_fields.get("jobid")
        if job_id:
            job = ensure_job(job_id)
            job.status = "failed"
            job.reason = error_fields.get("message", job.reason)
            if "external_jobid" in error_fields:
                job.slurm_jobid = error_fields["external_jobid"]
            if "log" in error_fields:
                job.log = _merge_log_values(job.log, error_fields["log"])
        message = error_fields.get("message") or error_fields.get("rule")
        if message:
            run.errors.append(message)
            last_event = message
        event_state = "failed"
        error_mode = False
        error_fields = {}

    for raw_line in path.read_text(errors="replace").splitlines():
        line = raw_line.rstrip()
        stripped = line.strip()
        if not stripped:
            finish_error_block()
            current_job = None
            continue
        if stripped.startswith(
            (
                "INFO:snakemake.logging:",
                "WARNING:snakemake.logging:",
                "Activating conda environment:",
            )
        ):
            continue
        if stripped.startswith("/user/") and "bind:" in stripped:
            continue

        if error_mode and not line.startswith("    "):
            finish_error_block()
        if current_job and not line.startswith("    "):
            current_job = None

        if match := _SLURM_RUN_ID_RE.match(stripped):
            run.slurm_run_id = match.group(1)
            last_event = stripped
            continue

        if match := _PROGRESS_RE.search(stripped):
            run.progress_done = int(match.group(1))
            run.progress_total = int(match.group(2))
            run.progress_percent = int(match.group(3))
            last_event = stripped
            event_state = "running"
            continue

        if match := _RULE_RE.match(stripped):
            current_job = {"kind": match.group(1), "rule": match.group(2)}
            continue

        if current_job and (match := _FIELD_RE.match(line)):
            key, value = match.groups()
            current_job[key.replace("-", "_")] = value
            if "jobid" in current_job:
                job = ensure_job(current_job["jobid"])
                job.kind = current_job.get("kind", job.kind)
                job.rule = current_job.get("rule", job.rule)
                job.wildcards = current_job.get("wildcards", job.wildcards)
                job.input = current_job.get("input", job.input)
                job.output = current_job.get("output", job.output)
                job.log = _merge_log_values(job.log, current_job.get("log", ""))
                job.benchmark = current_job.get("benchmark", job.benchmark)
                if job.kind == "localrule":
                    job.status = "running"
            continue

        if match := _SUBMIT_RE.search(stripped):
            snk_jobid, slurm_jobid, slurm_log = match.groups()
            job = ensure_job(snk_jobid)
            job.slurm_jobid = slurm_jobid
            job.slurm_log = slurm_log
            job.status = "submitted"
            last_event = f"Submitted {job.rule or snk_jobid} as {slurm_jobid}"
            event_state = "running"
            current_job = None
            continue

        if match := _FINISH_RE.search(stripped):
            snk_jobid, rule = match.groups()
            job = ensure_job(snk_jobid)
            job.rule = job.rule or rule
            job.status = "finished"
            if job.kind == "localrule" and not job.slurm_jobid:
                job.slurm_state = "LOCAL"
            last_event = f"Finished {job.rule or snk_jobid}"
            event_state = "running"
            current_job = None
            continue

        if match := _ERROR_RULE_RE.match(stripped):
            current_job = None
            error_mode = True
            error_fields = {"rule": match.group(1)}
            continue

        if error_mode and (match := _FIELD_RE.match(line)):
            key, value = match.groups()
            error_fields[key.replace("-", "_")] = value
            continue

        if stripped.startswith(
            (
                "LockException:",
                "WorkflowError:",
                "ProtectedOutputException",
                "MissingOutputException",
            )
        ):
            run.errors.append(stripped)
            last_event = stripped
            event_state = "failed"
            continue

        if stripped.startswith("Nothing to be done"):
            last_event = stripped
            event_state = "success"
            continue

        if stripped.startswith("Complete log(s):") or stripped == "Report created.":
            last_event = stripped
            event_state = "success"
            continue

    finish_error_block()
    run.state = event_state
    run.status_line = last_event
    return run


def apply_squeue_state(run: ControllerRun, squeue_output: str) -> None:
    entries: dict[str, tuple[str, str, str]] = {}
    for line in squeue_output.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split("|", maxsplit=4)
        if len(parts) < 5:
            continue
        job_id, state, elapsed, node, reason = parts
        entries[job_id] = (state, elapsed, node, reason)

    for job in run.jobs.values():
        if not job.slurm_jobid:
            continue
        if job.slurm_jobid not in entries:
            continue
        state, elapsed, node, reason = entries[job.slurm_jobid]
        job.slurm_state = state
        job.elapsed = elapsed
        job.node = node
        job.reason = reason if reason != "None" else job.reason
        if state in {"RUNNING", "COMPLETING"}:
            job.status = "running"
        elif state in {"PENDING", "CONFIGURING"}:
            job.status = "pending"
        else:
            job.status = state.lower()


def fetch_squeue_output(*, user: str, host: str = "") -> str:
    if not user:
        return ""

    command = f"squeue -h -u {user} -o '%i|%T|%M|%N|%R'"
    if host:
        with Connection(host=host, user=user) as connection:
            result = connection.run(command, hide=True, warn=True)
            return result.stdout if result.ok else ""

    result = subprocess.run(
        command,
        shell=True,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return ""
    return result.stdout


def _job_sort_key(job: JobRecord) -> tuple[int, str, str]:
    priority = {
        "running": 0,
        "pending": 1,
        "submitted": 2,
        "failed": 3,
        "finished": 4,
    }
    return (
        priority.get(job.status, 5),
        job.rule,
        job.snakemake_id.zfill(8),
    )


def _iter_display_jobs(run: ControllerRun, limit: int) -> Iterable[JobRecord]:
    jobs = sorted(run.jobs.values(), key=_job_sort_key)
    if limit <= 0:
        return jobs
    return jobs[:limit]


def job_rows(run: ControllerRun, *, limit: int = 0) -> list[dict[str, str]]:
    rows = []
    for job in _iter_display_jobs(run, limit):
        state = job.status
        if job.slurm_state:
            state = f"{state}/{job.slurm_state}"
        log_value = job.log or job.slurm_log or "-"
        display_log = log_value.split(",")[0].strip() if log_value != "-" else "-"
        display_log = re.sub(
            r"\s*\(check log file\(s\) for error details\)\s*$",
            "",
            display_log,
        )
        rows.append(
            {
                "smk_jobid": job.snakemake_id,
                "slurm_jobid": job.slurm_jobid or "-",
                "state": state,
                "rule": job.rule or "-",
                "wildcards": job.wildcards or "-",
                "elapsed": job.elapsed or "-",
                "node_or_reason": job.node or job.reason or "-",
                "input": job.input or "-",
                "output": job.output or "-",
                "log": log_value,
                "log_file": display_log,
            }
        )
    return rows


def resolve_log_paths(project_dir: Path, log_value: str) -> list[Path]:
    if not log_value or log_value == "-":
        return []

    paths: list[Path] = []
    seen: set[Path] = set()
    for raw_part in log_value.split(","):
        cleaned = raw_part.strip()
        cleaned = re.sub(r"\s*\(check log file\(s\) for error details\)\s*$", "", cleaned)
        if not cleaned:
            continue
        path = Path(cleaned)
        if not path.is_absolute():
            path = project_dir / path
        if path in seen:
            continue
        seen.add(path)
        paths.append(path)
    return paths


def _normalize_display_lines(lines: Iterable[str]) -> list[str]:
    normalized: list[str] = []
    previous_line: str | None = None
    for raw_line in lines:
        line = raw_line.rstrip()
        if line.startswith(_DISPLAY_LOG_PREFIXES):
            continue
        if line.startswith("/user/") and "bind:" in line:
            continue
        if previous_line is not None and line == previous_line:
            continue
        if not line and previous_line == "":
            continue
        normalized.append(line)
        previous_line = line
    return normalized


def format_log_text(
    text: str,
    *,
    source: str | Path,
    note: str = "",
    max_lines: int = 300,
    max_chars: int = 50000,
    normalize: bool = True,
) -> str:
    lines = text.splitlines()
    if normalize:
        lines = _normalize_display_lines(lines)

    tail_lines = lines[-max_lines:]
    tail = "\n".join(tail_lines)
    if len(tail) > max_chars:
        tail = tail[-max_chars:]

    header_lines = [f"# {source}"]
    if note:
        header_lines.append(f"# {note}")
    if len(lines) > max_lines:
        header_lines.append(f"# Showing last {max_lines} of {len(lines)} lines.")

    header = "\n".join(header_lines) + "\n\n"
    return header + tail


def _block_field_value(block_lines: Iterable[str], field_name: str) -> str:
    normalized_field = field_name.replace("-", "_")
    for line in block_lines:
        match = _FIELD_RE.match(line)
        if not match:
            continue
        key, value = match.groups()
        if key.replace("-", "_") == normalized_field:
            return value
    return ""


def _collect_indented_block(lines: list[str], start: int) -> tuple[list[str], int]:
    block = [lines[start]]
    index = start + 1
    while index < len(lines):
        line = lines[index]
        if line.startswith("    "):
            block.append(line)
            index += 1
            continue
        if not line.strip():
            block.append(line)
            index += 1
            continue
        break
    return block, index


def controller_job_log_context(
    controller_log: Path,
    job: JobRecord,
    *,
    max_sections: int = 4,
    max_lines: int = 200,
    max_chars: int = 30000,
) -> list[tuple[str, str]]:
    if not controller_log.exists():
        return []

    try:
        lines = controller_log.read_text(errors="replace").splitlines()
    except Exception as exc:  # pragma: no cover - best effort viewer
        return [
            (
                "Controller context",
                format_log_text(
                    f"Unable to read file: {exc}",
                    source=controller_log,
                    note=f"Filtered controller context for Snakemake job {job.snakemake_id}",
                    max_lines=max_lines,
                    max_chars=max_chars,
                    normalize=False,
                ),
            )
        ]

    raw_sections: list[tuple[str, list[str]]] = []
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()

        if _RULE_RE.match(stripped):
            block, next_index = _collect_indented_block(lines, index)
            if _block_field_value(block, "jobid") == job.snakemake_id:
                section = list(block)
                index = next_index
                while index < len(lines):
                    current_line = lines[index]
                    current_stripped = current_line.strip()
                    if _RULE_RE.match(current_stripped) or _ERROR_RULE_RE.match(current_stripped):
                        break
                    submit_match = _SUBMIT_RE.search(current_stripped)
                    finish_match = _FINISH_RE.search(current_stripped)
                    section.append(current_line)
                    index += 1
                    if submit_match and submit_match.group(1) == job.snakemake_id:
                        break
                    if finish_match and finish_match.group(1) == job.snakemake_id:
                        break
                raw_sections.append(("Controller events", section))
                continue
            index = next_index
            continue

        if _ERROR_RULE_RE.match(stripped):
            block, next_index = _collect_indented_block(lines, index)
            error_jobid = _block_field_value(block, "jobid")
            external_jobid = _block_field_value(block, "external_jobid")
            if error_jobid == job.snakemake_id or (
                job.slurm_jobid and external_jobid == job.slurm_jobid
            ):
                raw_sections.append(("Controller error", block))
            index = next_index
            continue

        index += 1

    sections: list[tuple[str, str]] = []
    seen_contents: set[str] = set()
    for title, section_lines in raw_sections:
        normalized_text = "\n".join(_normalize_display_lines(section_lines)).strip()
        if not normalized_text or normalized_text in seen_contents:
            continue
        seen_contents.add(normalized_text)
        sections.append(
            (
                title,
                format_log_text(
                    normalized_text,
                    source=controller_log,
                    note=f"{title} for Snakemake job {job.snakemake_id}",
                    max_lines=max_lines,
                    max_chars=max_chars,
                    normalize=False,
                ),
            )
        )
        if len(sections) >= max_sections:
            break
    return sections


def read_log_tail(path: Path, *, max_lines: int = 300, max_chars: int = 50000) -> str:
    if not path.exists():
        return f"{path}\n\nFile not found."

    try:
        text = path.read_text(errors="replace")
    except Exception as exc:  # pragma: no cover - best effort viewer
        return f"{path}\n\nUnable to read file: {exc}"

    return format_log_text(text, source=path, max_lines=max_lines, max_chars=max_chars)


def render_monitor(run: ControllerRun, *, workflow_name: str, limit: int = 20) -> None:
    console = Console()
    counts = run.counts()
    progress = "-"
    if run.progress_total:
        progress = f"{run.progress_done}/{run.progress_total} ({run.progress_percent}%)"

    header = Table.grid(expand=False)
    header.add_column(style="bold cyan")
    header.add_column()
    header.add_row("Workflow", workflow_name)
    header.add_row("Controller", str(run.controller_log.relative_to(run.controller_log.parent.parent.parent)))
    header.add_row("Started", run.started_at)
    header.add_row("State", run.state)
    header.add_row("Progress", progress)
    if run.slurm_run_id:
        header.add_row("SLURM Run", run.slurm_run_id)
    header.add_row("Last Event", run.status_line)
    header.add_row(
        "Counts",
        (
            f"running={counts['running']}  pending={counts['submitted'] - counts['running']}  "
            f"finished={counts['finished']}  failed={counts['failed']}  local={counts['local']}"
        ),
    )
    console.print(header)

    if run.errors:
        console.print()
        console.print("[bold red]Controller Errors[/bold red]")
        for message in run.errors[-3:]:
            console.print(f"- {message}")

    jobs = list(_iter_display_jobs(run, limit))
    if not jobs:
        return

    table = Table(title="Jobs", show_lines=False)
    table.add_column("SMK", justify="right")
    table.add_column("SLURM", justify="right")
    table.add_column("State")
    table.add_column("Rule")
    table.add_column("Wildcards", overflow="fold")
    table.add_column("Elapsed")
    table.add_column("Node/Reason", overflow="fold")

    for job in jobs:
        state = job.status
        if job.slurm_state:
            state = f"{state}/{job.slurm_state}"
        table.add_row(
            job.snakemake_id,
            job.slurm_jobid or "-",
            state,
            job.rule or "-",
            job.wildcards or "-",
            job.elapsed or "-",
            job.node or job.reason or "-",
        )
    console.print()
    console.print(table)
