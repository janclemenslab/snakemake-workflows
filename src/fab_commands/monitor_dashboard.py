"""FastHTML dashboard for Snakemake controller logs."""

from __future__ import annotations

import os
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fasthtml.common import *  # noqa: F403
from monsterui.all import *  # noqa: F403

from fab_commands.monitor import (
    apply_squeue_state,
    controller_job_log_context,
    controller_logs,
    fetch_squeue_output,
    job_rows,
    parse_controller_log,
    read_log_tail,
    resolve_log_paths,
)


PROJECT_DIR = Path(os.environ.get("FAB_MONITOR_PROJECT_DIR", ".")).resolve()
DEFAULT_USER = os.environ.get("FAB_MONITOR_USER", "").strip()
DEFAULT_REMOTE_HOST = os.environ.get("FAB_MONITOR_REMOTE_HOST", "").strip()
SORTABLE_COLUMNS = {
    "state": "state",
    "rule": "rule",
    "wildcards": "wildcards",
    "elapsed": "elapsed",
    "node_or_reason": "node_or_reason",
    "log_file": "log_file",
}


def _env_int(name: str, default: int, *, minimum: int = 0) -> int:
    try:
        value = int(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = default
    return max(minimum, value)


DEFAULT_LIMIT = _env_int("FAB_MONITOR_LIMIT", 100, minimum=10)
DEFAULT_REFRESH = _env_int("FAB_MONITOR_REFRESH", 5, minimum=0)
DEFAULT_HOST = os.environ.get("FAB_MONITOR_BIND_HOST", "127.0.0.1").strip() or "127.0.0.1"
DEFAULT_PORT = _env_int("FAB_MONITOR_PORT", 2718, minimum=1)
PAGE_SIZE_OPTIONS = [50, 100, 200]


app, rt = fast_app(  # noqa: F405
    pico=False,
    live=False,
    title="Snakemake Monitor",
    hdrs=Theme.slate.headers(),  # noqa: F405
)


def _available_logs() -> list[Path]:
    return controller_logs(PROJECT_DIR)


def _normalize_limit(limit: int | str | None) -> int:
    try:
        value = int(limit if limit is not None else DEFAULT_LIMIT)
    except (TypeError, ValueError):
        value = DEFAULT_LIMIT
    return max(10, min(value, 200))


def _normalize_page(page: int | str | None) -> int:
    try:
        value = int(page if page is not None else 1)
    except (TypeError, ValueError):
        value = 1
    return max(1, value)


def _normalize_refresh(refresh: int | str | None) -> int:
    try:
        value = int(refresh if refresh is not None else DEFAULT_REFRESH)
    except (TypeError, ValueError):
        value = DEFAULT_REFRESH
    return max(0, min(value, 60))


def _selected_log_path(controller_log: str | None) -> str:
    logs = _available_logs()
    if controller_log:
        candidate = Path(controller_log).expanduser()
        if candidate.exists():
            return str(candidate.resolve())
    if logs:
        return str(logs[-1])
    return ""


def _load_run(controller_log: str, user: str, remote_host: str):
    selected_log = _selected_log_path(controller_log)
    if not selected_log:
        return None, "", ""

    run = parse_controller_log(Path(selected_log))
    queue_warning = ""
    queue_user = (user or "").strip()
    queue_host = (remote_host or "").strip()

    if queue_user and queue_host:
        squeue_output = fetch_squeue_output(user=queue_user, host=queue_host)
        if squeue_output:
            apply_squeue_state(run, squeue_output)
        else:
            queue_warning = "Could not fetch live SLURM state."
    elif queue_user and not queue_host:
        queue_warning = "Set an SSH host to enrich jobs with live SLURM state."

    return run, queue_warning, selected_log


def _elapsed_sort_value(value: str) -> int:
    if not value or value == "-":
        return -1
    days = 0
    clock = value
    if "-" in value:
        day_text, clock = value.split("-", maxsplit=1)
        if day_text.isdigit():
            days = int(day_text)
    parts = clock.split(":")
    if len(parts) != 3 or not all(part.isdigit() for part in parts):
        return -1
    hours, minutes, seconds = (int(part) for part in parts)
    return (((days * 24) + hours) * 60 + minutes) * 60 + seconds


def _sort_value(row: dict[str, str], sort_by: str):
    if sort_by == "elapsed":
        return (_elapsed_sort_value(row["elapsed"]), row["smk_jobid"])
    return (row.get(sort_by, "").lower(), row["smk_jobid"])


def _sorted_rows(run, sort_by: str, sort_dir: str):
    rows = job_rows(run, limit=0)
    if sort_by in SORTABLE_COLUMNS:
        rows = sorted(rows, key=lambda row: _sort_value(row, sort_by), reverse=sort_dir == "desc")
    return rows


def _warning_card(title: str, messages: list[str], tone: str):
    tone_cls = {
        "warn": "border-orange-300 bg-orange-50",
        "danger": "border-red-300 bg-red-50",
    }.get(tone, "border-slate-300 bg-slate-50")
    return Card(  # noqa: F405
        Ul(*[Li(message) for message in messages], cls="list-disc pl-5 space-y-1"),  # noqa: F405
        header=H4(title, cls="text-sm font-semibold"),  # noqa: F405
        cls=f"h-full border {tone_cls}",
        body_cls="space-y-3",
    )


def _controls(
    controller_log: str,
    user: str,
    remote_host: str,
    limit: int,
    page: int,
    refresh: int,
    sort_by: str,
    sort_dir: str,
    controls_open: bool = False,
):
    logs = list(reversed(_available_logs()))
    run_options = [
        Option(path.name, value=str(path), selected=str(path) == controller_log)  # noqa: F405
        for path in logs
    ]

    refresh_options = [
        ("Off", 0),
        ("3s", 3),
        ("5s", 5),
        ("10s", 10),
        ("30s", 30),
    ]

    run_select = Div(
        DivFullySpaced(  # noqa: F405
            P("Run", cls="text-sm font-medium"),
            Button(  # noqa: F405
                UkIcon("refresh-cw", width=14, height=14),  # noqa: F405
                submit=False,
                cls=(ButtonT.ghost, ButtonT.sm),  # noqa: F405
                title="Refresh controller log list",
                hx_get="/controls?controls_open=1",
                hx_target="#controls-card",
                hx_swap="outerHTML",
                hx_include="#controls-form",
            ),
        ),
        Select(  # noqa: F405
            *run_options,
            id="controller_log",
            name="controller_log",
            placeholder="No controller logs",
            searchable=False,
            disabled=not run_options,
        ),
        cls="space-y-2",
    )

    refresh_select = Div(
        P("Refresh", cls="text-sm font-medium"),
        Select(  # noqa: F405
            *[
                Option(label, value=str(seconds), selected=seconds == refresh)  # noqa: F405
                for label, seconds in refresh_options
            ],
            id="refresh",
            name="refresh",
            searchable=False,
        ),
        cls="space-y-2",
    )

    form = Form(  # noqa: F405
        Input(type="hidden", id="page", name="page", value=str(page)),  # noqa: F405
        Input(type="hidden", id="sort_by", name="sort_by", value=sort_by),  # noqa: F405
        Input(type="hidden", id="sort_dir", name="sort_dir", value=sort_dir),  # noqa: F405
        Input(type="hidden", id="limit", name="limit", value=str(limit)),  # noqa: F405
        Grid(  # noqa: F405
            run_select,
            LabelInput("User", id="user", value=user, placeholder="abcd1234"),  # noqa: F405
            LabelInput(  # noqa: F405
                "Host",
                id="remote_host",
                value=remote_host,
                placeholder="rosa.hpc.uni-oldenburg.de",
            ),
            refresh_select,
            Div(
                Button(  # noqa: F405
                    "Apply",
                    submit=False,
                    cls=(ButtonT.primary, ButtonT.sm),  # noqa: F405
                    onclick="document.getElementById('page').value='1';",
                    hx_get="/dashboard-shell",
                    hx_target="#dashboard-shell",
                    hx_swap="outerHTML",
                    hx_include="#controls-form",
                ),
                cls="flex items-end h-full",
            ),
            cols_min=1,
            cols_md=2,
            cols_lg=3,
            cols_xl=5,
            cls="gap-2 items-end pt-2",
        ),
        id="controls-form",
        cls="space-y-0",
    )

    summary_line = DivFullySpaced(  # noqa: F405
        Div(
            H2(PROJECT_DIR.name, cls="text-lg font-semibold tracking-tight"),  # noqa: F405
            P(Path(controller_log).name if controller_log else "No controller log selected", cls="text-xs opacity-70"),
            cls="space-y-0.5",
        ),
        P("Options", cls="text-sm opacity-60"),
    )

    return Card(  # noqa: F405
        Details(  # noqa: F405
            Summary(summary_line, cls="cursor-pointer list-none"),  # noqa: F405
            form,
            open=controls_open,
            cls="group",
        ),
        id="controls-card",
        body_cls="space-y-0 py-3",
    )


def _summary(run, queue_warning: str, selected_log: str):
    if run is None:
        return Card(  # noqa: F405
            P("No controller logs found under `log/slurm`.", cls="text-sm"),
            cls="border border-amber-300 bg-amber-50",
            body_cls="space-y-0 py-3",
        )

    counts = run.counts()
    progress = "-"
    if run.progress_total:
        progress = f"{run.progress_done}/{run.progress_total} ({run.progress_percent}%)"

    compact_row = DivFullySpaced(  # noqa: F405
        Div(
            H3(Path(selected_log).name, cls="text-base font-semibold"),  # noqa: F405
            P(
                f"{run.started_at} | {run.state} | {progress}",
                cls="text-sm opacity-75",
            ),
            cls="space-y-0.5",
        ),
        Div(
            Span("running ", Strong(str(counts["running"]))),  # noqa: F405
            Span("pending ", Strong(str(counts["submitted"] - counts["running"]))),  # noqa: F405
            Span("failed ", Strong(str(counts["failed"]))),  # noqa: F405
            cls="flex flex-wrap gap-4 text-sm whitespace-nowrap",
        ),
    )

    blocks = [compact_row]
    meta = []
    if run.slurm_run_id:
        meta.append(f"SLURM run ID: {run.slurm_run_id}")
    if run.status_line:
        meta.append(run.status_line)
    if meta:
        blocks.append(P(" | ".join(meta), cls="text-sm opacity-70 break-all"))
    warnings = []
    if queue_warning:
        warnings.append(_warning_card("Live Queue", [queue_warning], "warn"))
    error_messages = [message for message in run.errors if message and message != "None"]
    if error_messages:
        warnings.append(_warning_card("Recent Errors", error_messages[-5:], "danger"))
    if warnings:
        blocks.append(Grid(*warnings, cols_min=1, cols_lg=2, cls="gap-2"))  # noqa: F405
    return Card(*blocks, body_cls="space-y-2 py-3")  # noqa: F405


def _log_link(row: dict[str, str]):
    if row["log_file"] == "-":
        return CodeSpan("-", cls="text-xs")  # noqa: F405
    return A(  # noqa: F405
        row["log_file"],
        href="#",
        cls="font-mono text-xs hover:underline",
        title=row["log_file"],
        hx_get=f"/job-modal?smk_jobid={row['smk_jobid']}",
        hx_target="#modal-container",
        hx_swap="innerHTML",
        hx_include="#controls-form",
    )


def _sort_header(label: str, key: str, sort_by: str, sort_dir: str):
    marker = ""
    next_dir = "asc"
    if sort_by == key:
        marker = " ^" if sort_dir == "asc" else " v"
        next_dir = "desc" if sort_dir == "asc" else "asc"
    return Th(  # noqa: F405
        A(  # noqa: F405
            f"{label}{marker}",
            href="#",
            cls="hover:underline whitespace-nowrap",
            onclick=(
                f"document.getElementById('sort_by').value='{key}';"
                f"document.getElementById('sort_dir').value='{next_dir}';"
            ),
            hx_get="/dashboard-shell",
            hx_target="#dashboard-shell",
            hx_swap="outerHTML",
            hx_include="#controls-form",
        )
    )


def _pager_button(label: str, page: int, disabled: bool = False):
    classes = (ButtonT.secondary, ButtonT.sm)
    if disabled:
        return Button(label, submit=False, cls=classes, disabled=True)  # noqa: F405
    return Button(  # noqa: F405
        label,
        submit=False,
        cls=classes,
        onclick=f"document.getElementById('page').value='{page}';",
        hx_get="/dashboard-shell",
        hx_target="#dashboard-shell",
        hx_swap="outerHTML",
        hx_include="#controls-form",
    )


def _page_size_button(size: int, current: int):
    classes = (ButtonT.primary if size == current else ButtonT.ghost, ButtonT.sm)
    return Button(  # noqa: F405
        str(size),
        submit=False,
        cls=classes,
        onclick=(
            f"document.getElementById('limit').value='{size}';"
            "document.getElementById('page').value='1';"
        ),
        hx_get="/dashboard-shell",
        hx_target="#dashboard-shell",
        hx_swap="outerHTML",
        hx_include="#controls-form",
    )


def _jobs_table(run, limit: int, page: int, sort_by: str, sort_dir: str):
    if run is None:
        return Div()  # noqa: F405

    rows = _sorted_rows(run, sort_by, sort_dir)
    if not rows:
        return Card(P("No jobs were parsed from the selected controller log.", cls="text-sm"))  # noqa: F405

    total_rows = len(rows)
    page_count = max(1, (total_rows + limit - 1) // limit)
    current_page = min(max(1, page), page_count)
    start = (current_page - 1) * limit
    stop = start + limit
    visible_rows = rows[start:stop]

    body = []
    for row in visible_rows:
        body.append(
            Tr(  # noqa: F405
                Td(CodeSpan(row["state"]), cls="whitespace-nowrap text-xs"),  # noqa: F405
                Td(row["rule"], cls="font-medium"),
                Td(row["wildcards"], cls="max-w-xl text-sm"),
                Td(row["elapsed"], cls="whitespace-nowrap text-sm"),
                Td(row["node_or_reason"], cls="max-w-sm text-sm"),
                Td(_log_link(row), cls="max-w-sm"),
            )
        )

    return Card(  # noqa: F405
        Div(
            Table(  # noqa: F405
                Thead(  # noqa: F405
                    Tr(  # noqa: F405
                        _sort_header("State", "state", sort_by, sort_dir),
                        _sort_header("Rule", "rule", sort_by, sort_dir),
                        _sort_header("Wildcards", "wildcards", sort_by, sort_dir),
                        _sort_header("Elapsed", "elapsed", sort_by, sort_dir),
                        _sort_header("Node / Reason", "node_or_reason", sort_by, sort_dir),
                        _sort_header("Log", "log_file", sort_by, sort_dir),
                    )
                ),
                Tbody(*body),  # noqa: F405
            ),
            cls="overflow-x-auto",
        ),
        footer=DivFullySpaced(  # noqa: F405
            Div(
                _pager_button("Previous", current_page - 1, disabled=current_page <= 1),
                P(
                    f"Page {current_page} of {page_count} | {total_rows} jobs",
                    cls="text-xs opacity-70",
                ),
                _pager_button("Next", current_page + 1, disabled=current_page >= page_count),
                cls="flex items-center gap-2 flex-wrap",
            ),
            Div(
                P("Rows per page", cls="text-xs opacity-70"),
                *[_page_size_button(size, limit) for size in PAGE_SIZE_OPTIONS],
                cls="flex items-center gap-2 flex-wrap",
            ),
        ),
        header=DivFullySpaced(  # noqa: F405
            H3("Jobs", cls="text-base font-semibold"),  # noqa: F405
            P(
                f"Showing {start + 1}-{min(stop, total_rows)} of {total_rows}",
                cls="text-xs opacity-70",
            ),
        ),
        body_cls="space-y-0 py-2",
    )


def _log_card(title: str, subtitle: str, content: str, *, cls: str = ""):
    return Card(  # noqa: F405
        Pre(  # noqa: F405
            content,
            cls="w-full overflow-auto rounded border p-3",
            style=(
                "height: 82vh; min-height: 82vh; overflow: auto; "
                "white-space: pre-wrap; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; "
                "font-size: 12px; line-height: 1.5;"
            ),
        ),
        header=Div(
            P(title, cls="text-sm font-medium"),
            P(subtitle, cls="text-xs opacity-70 break-all"),
            cls="space-y-1",
        ),
        body_cls="space-y-3",
        cls=cls,
    )


def _dashboard_shell(
    controller_log: str,
    user: str,
    remote_host: str,
    limit: int,
    page: int,
    refresh: int,
    sort_by: str,
    sort_dir: str,
):
    run, queue_warning, selected_log = _load_run(controller_log, user, remote_host)
    poll_attrs = {}
    if refresh > 0:
        poll_attrs = {
            "hx_get": "/dashboard-shell",
            "hx_trigger": f"load, every {refresh}s",
            "hx_target": "this",
            "hx_swap": "outerHTML",
            "hx_include": "#controls-form",
        }

    return Div(  # noqa: F405
        _summary(run, queue_warning, selected_log),
        _jobs_table(run, limit, page, sort_by, sort_dir),
        id="dashboard-shell",
        cls="space-y-2",
        **poll_attrs,
    )


@rt("/")  # noqa: F405
def index(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    page: int = 1,
    refresh: int = DEFAULT_REFRESH,
    sort_by: str = "",
    sort_dir: str = "asc",
):
    limit = _normalize_limit(limit)
    page = _normalize_page(page)
    refresh = _normalize_refresh(refresh)
    controller_log = _selected_log_path(controller_log)

    return Title(f"Snakemake Monitor | {PROJECT_DIR.name}"), Container(  # noqa: F405
        _controls(controller_log, user, remote_host, limit, page, refresh, sort_by, sort_dir),
        _dashboard_shell(controller_log, user, remote_host, limit, page, refresh, sort_by, sort_dir),
        Div(id="modal-container"),  # noqa: F405
        cls="space-y-3 py-3",
    )


@rt("/controls")  # noqa: F405
def controls(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    page: int = 1,
    refresh: int = DEFAULT_REFRESH,
    sort_by: str = "",
    sort_dir: str = "asc",
    controls_open: int = 0,
):
    return _controls(
        _selected_log_path(controller_log),
        user,
        remote_host,
        _normalize_limit(limit),
        _normalize_page(page),
        _normalize_refresh(refresh),
        sort_by,
        sort_dir,
        bool(int(controls_open)),
    )


@rt("/dashboard-shell")  # noqa: F405
def dashboard_shell(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    page: int = 1,
    refresh: int = DEFAULT_REFRESH,
    sort_by: str = "",
    sort_dir: str = "asc",
):
    return _dashboard_shell(
        _selected_log_path(controller_log),
        user,
        remote_host,
        _normalize_limit(limit),
        _normalize_page(page),
        _normalize_refresh(refresh),
        sort_by,
        sort_dir,
    )


@rt("/job-modal")  # noqa: F405
def job_modal(
    smk_jobid: str,
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
):
    run, _, _ = _load_run(controller_log, user, remote_host)
    if run is None:
        return Div()  # noqa: F405

    row_lookup = {row["smk_jobid"]: row for row in job_rows(run, limit=0)}
    row = row_lookup.get(smk_jobid)
    job = run.jobs.get(smk_jobid)
    if row is None or job is None:
        return Div()  # noqa: F405

    log_cards = []
    for title, content in controller_job_log_context(
        run.controller_log,
        job,
        max_lines=200,
        max_chars=40000,
    ):
        log_cards.append(
            _log_card(
                title,
                str(run.controller_log),
                content,
                cls="border border-sky-200 bg-sky-50",
            )
        )

    for path in resolve_log_paths(PROJECT_DIR, row["log"]):
        log_cards.append(
            _log_card(
                path.name,
                str(path),
                read_log_tail(path, max_lines=250, max_chars=40000),
            )
        )

    if not log_cards:
        log_cards = [
            Card(
                P("No readable log files for this job.", cls="text-sm"),
                cls="border border-amber-300 bg-amber-50",
            )
        ]

    header = Div(
        H3(row["rule"], cls="text-lg font-semibold"),  # noqa: F405
        Div(
            CodeSpan(row["state"]),  # noqa: F405
            P(f"SMK {row['smk_jobid']} | SLURM {row['slurm_jobid']}", cls="text-sm opacity-70"),
            cls="space-y-1",
        ),
        P(row["wildcards"], cls="text-sm break-all"),
        cls="space-y-2",
    )

    return Modal(  # noqa: F405
        *log_cards,
        header=header,
        footer=Div(
            ModalCloseButton("Close", cls=(ButtonT.secondary, "static"), submit=False),  # noqa: F405
            cls="flex justify-end",
        ),
        id="job-log-modal",
        hx_init=True,
        hx_open=True,
        dialog_cls="uk-width-5-6@xl uk-width-11-12 max-h-[92vh]",
        body_cls="space-y-4 overflow-y-auto",
        style="padding-top: 2vh; padding-bottom: 2vh;",
    )


serve(host=DEFAULT_HOST, port=DEFAULT_PORT, reload=False)  # noqa: F405
