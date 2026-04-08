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


def _env_int(name: str, default: int, *, minimum: int = 0) -> int:
    try:
        value = int(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = default
    return max(minimum, value)


DEFAULT_LIMIT = _env_int("FAB_MONITOR_LIMIT", 25, minimum=5)
DEFAULT_REFRESH = _env_int("FAB_MONITOR_REFRESH", 5, minimum=0)
DEFAULT_HOST = os.environ.get("FAB_MONITOR_BIND_HOST", "127.0.0.1").strip() or "127.0.0.1"
DEFAULT_PORT = _env_int("FAB_MONITOR_PORT", 2718, minimum=1)


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
    return max(5, min(value, 100))


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


def _stats_card(label: str, value: str | int, caption: str = ""):
    body = [Div(str(value), cls="text-2xl font-semibold leading-none")]
    if caption:
        body.append(P(caption, cls="mt-2 text-sm opacity-70"))
    return Card(  # noqa: F405
        *body,
        header=P(label, cls="text-xs font-medium uppercase tracking-wide opacity-70"),
        body_cls="space-y-0",
        cls="h-full",
    )


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


def _controls(controller_log: str, user: str, remote_host: str, limit: int, refresh: int):
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
        P("Run", cls="text-sm font-medium"),
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

    return Card(  # noqa: F405
        Form(  # noqa: F405
            Grid(  # noqa: F405
                run_select,
                LabelInput("User", id="user", value=user, placeholder="abcd1234"),  # noqa: F405
                LabelInput(  # noqa: F405
                    "Host",
                    id="remote_host",
                    value=remote_host,
                    placeholder="rosa.hpc.uni-oldenburg.de",
                ),
                LabelInput("Rows", id="limit", value=str(limit), type="number", min="5", max="100"),  # noqa: F405
                refresh_select,
                Div(
                    Button(  # noqa: F405
                        "Apply",
                        submit=False,
                        cls=ButtonT.primary,  # noqa: F405
                        hx_get="/dashboard-shell",
                        hx_target="#dashboard-shell",
                        hx_swap="outerHTML",
                        hx_include="#controls-form",
                    ),
                    cls="flex items-end h-full",
                ),
                cols=1,
                cols_md=2,
                cols_lg=3,
                cols_xl=6,
                cls="gap-3 items-end",
            ),
            id="controls-form",
            cls="space-y-0",
        ),
        header=DivFullySpaced(  # noqa: F405
            Div(
                H2(PROJECT_DIR.name, cls="text-xl font-semibold tracking-tight"),  # noqa: F405
                P("Snakemake controller monitor", cls="text-sm opacity-70"),
                cls="space-y-1",
            ),
            P("FastHTML dashboard", cls="text-sm opacity-60"),
        ),
        body_cls="space-y-0",
    )


def _summary(run, queue_warning: str):
    if run is None:
        return Card(  # noqa: F405
            P("No controller logs found under `log/slurm`.", cls="text-sm"),
            cls="border border-amber-300 bg-amber-50",
        )

    counts = run.counts()
    progress = "-"
    if run.progress_total:
        progress = f"{run.progress_done}/{run.progress_total} ({run.progress_percent}%)"

    items = [
        _stats_card("State", run.state, run.started_at),
        _stats_card("Progress", progress, run.status_line),
        _stats_card("Running", counts["running"]),
        _stats_card("Pending", counts["submitted"] - counts["running"]),
        _stats_card("Failed", counts["failed"]),
    ]

    blocks = [Grid(*items, cols=1, cols_md=2, cols_lg=3, cols_xl=5, cls="gap-3")]  # noqa: F405
    warnings = []
    if queue_warning:
        warnings.append(_warning_card("Live Queue", [queue_warning], "warn"))
    if run.errors:
        warnings.append(_warning_card("Recent Errors", run.errors[-5:], "danger"))
    if warnings:
        blocks.append(Grid(*warnings, cols=1, cols_lg=2, cls="gap-3"))  # noqa: F405
    return Div(*blocks, cls="space-y-3")  # noqa: F405


def _job_button(row: dict[str, str]):
    return Button(  # noqa: F405
        "View logs",
        submit=False,
        cls=ButtonT.secondary,  # noqa: F405
        hx_get=f"/job-modal?smk_jobid={row['smk_jobid']}",
        hx_target="#modal-container",
        hx_swap="innerHTML",
        hx_include="#controls-form",
    )


def _jobs_table(run, limit: int):
    if run is None:
        return Div()  # noqa: F405

    rows = job_rows(run, limit=limit)
    if not rows:
        return Card(P("No jobs were parsed from the selected controller log.", cls="text-sm"))  # noqa: F405

    body = []
    for row in rows:
        body.append(
            Tr(  # noqa: F405
                Td(_job_button(row), cls="whitespace-nowrap"),
                Td(CodeSpan(row["state"]), cls="whitespace-nowrap"),  # noqa: F405
                Td(row["rule"], cls="font-medium"),
                Td(row["wildcards"], cls="max-w-xl text-sm"),
                Td(row["elapsed"], cls="whitespace-nowrap text-sm"),
                Td(row["node_or_reason"], cls="max-w-sm text-sm"),
                Td(CodeSpan(row["log_file"]), cls="max-w-sm text-xs"),
            )
        )

    return Card(  # noqa: F405
        Div(
            Table(  # noqa: F405
                Thead(  # noqa: F405
                    Tr(  # noqa: F405
                        Th(""),  # noqa: F405
                        Th("State"),  # noqa: F405
                        Th("Rule"),  # noqa: F405
                        Th("Wildcards"),  # noqa: F405
                        Th("Elapsed"),  # noqa: F405
                        Th("Node / Reason"),  # noqa: F405
                        Th("Log"),  # noqa: F405
                    )
                ),
                Tbody(*body),  # noqa: F405
            ),
            cls="overflow-x-auto",
        ),
        header=DivFullySpaced(  # noqa: F405
            H3("Jobs", cls="text-lg font-semibold"),  # noqa: F405
            P(f"Showing up to {limit} jobs", cls="text-sm opacity-70"),
        ),
        body_cls="space-y-0",
    )


def _dashboard_shell(controller_log: str, user: str, remote_host: str, limit: int, refresh: int):
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

    headline = [
        P("Latest run" if selected_log else "Run", cls="text-xs uppercase tracking-wide opacity-60"),
        CodeSpan(Path(selected_log).name if selected_log else "No controller log selected"),  # noqa: F405
    ]
    if run and run.slurm_run_id:
        headline.append(P(f"SLURM run ID: {run.slurm_run_id}", cls="text-sm opacity-70"))

    return Div(  # noqa: F405
        Card(  # noqa: F405
            Div(*headline, cls="space-y-1"),
            body_cls="space-y-1",
        ),
        _summary(run, queue_warning),
        _jobs_table(run, limit),
        id="dashboard-shell",
        cls="space-y-3",
        **poll_attrs,
    )


@rt("/")  # noqa: F405
def index(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    refresh: int = DEFAULT_REFRESH,
):
    limit = _normalize_limit(limit)
    refresh = _normalize_refresh(refresh)
    controller_log = _selected_log_path(controller_log)

    return Title(f"Snakemake Monitor | {PROJECT_DIR.name}"), Container(  # noqa: F405
        _controls(controller_log, user, remote_host, limit, refresh),
        _dashboard_shell(controller_log, user, remote_host, limit, refresh),
        Div(id="modal-container"),  # noqa: F405
        cls="space-y-4 py-4",
    )


@rt("/dashboard-shell")  # noqa: F405
def dashboard_shell(
    controller_log: str = "",
    user: str = DEFAULT_USER,
    remote_host: str = DEFAULT_REMOTE_HOST,
    limit: int = DEFAULT_LIMIT,
    refresh: int = DEFAULT_REFRESH,
):
    return _dashboard_shell(
        _selected_log_path(controller_log),
        user,
        remote_host,
        _normalize_limit(limit),
        _normalize_refresh(refresh),
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
    if row is None:
        return Div()  # noqa: F405

    log_cards = []
    for path in resolve_log_paths(PROJECT_DIR, row["log"]):
        log_cards.append(
            Card(  # noqa: F405
                Textarea(  # noqa: F405
                    read_log_tail(path, max_lines=250, max_chars=40000),
                    readonly=True,
                    cls="min-h-[24rem] w-full font-mono text-xs leading-5",
                ),
                header=Div(
                    P(path.name, cls="text-sm font-medium"),
                    P(str(path), cls="text-xs opacity-70 break-all"),
                    cls="space-y-1",
                ),
                body_cls="space-y-3",
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
            ModalCloseButton("Close", cls=(ButtonT.secondary, "static")),  # noqa: F405
            cls="flex justify-end",
        ),
        id="job-log-modal",
        hx_init=True,
        hx_open=True,
        dialog_cls="uk-width-5-6@xl uk-width-11-12",
        body_cls="space-y-4",
    )


serve(host=DEFAULT_HOST, port=DEFAULT_PORT, reload=False)  # noqa: F405
