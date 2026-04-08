"""marimo dashboard for Snakemake controller logs."""

import marimo

__generated_with = "0.22.5"
app = marimo.App(width="full")


@app.cell
def _():
    import argparse
    import os
    from pathlib import Path

    import marimo as mo

    from fab_commands.monitor import (
        apply_squeue_state,
        controller_logs,
        fetch_squeue_output,
        job_rows,
        parse_controller_log,
        read_log_tail,
        resolve_log_paths,
    )

    return (
        Path,
        apply_squeue_state,
        argparse,
        controller_logs,
        fetch_squeue_output,
        job_rows,
        mo,
        os,
        parse_controller_log,
        read_log_tail,
        resolve_log_paths,
    )


@app.cell
def _(Path, argparse, os):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--user", default="")
    parser.add_argument("--remote-host", default="")
    parser.add_argument("--limit", type=int, default=25)
    args, _ = parser.parse_known_args()

    env_project_dir = os.environ.get("FAB_MONITOR_PROJECT_DIR", "").strip()
    env_user = os.environ.get("FAB_MONITOR_USER", "").strip()
    env_remote_host = os.environ.get("FAB_MONITOR_REMOTE_HOST", "").strip()
    env_limit = os.environ.get("FAB_MONITOR_LIMIT", "").strip()

    project_dir = Path(env_project_dir or args.project_dir).resolve()
    default_host = env_remote_host or args.remote_host.strip()
    default_user = env_user or args.user.strip()
    parsed_limit = int(env_limit) if env_limit.isdigit() else args.limit
    default_limit = max(5, parsed_limit)
    return default_host, default_limit, default_user, project_dir


@app.cell
def _(controller_logs, project_dir):
    logs = controller_logs(project_dir)
    options = {path.name: str(path) for path in reversed(logs)}
    default_log = logs[-1].name if logs else None
    return default_log, logs, options


@app.cell
def _(default_host, default_limit, default_log, default_user, mo, options, project_dir):
    controller_log = mo.ui.dropdown(
        options,
        value=default_log,
        label="Run",
    )
    user = mo.ui.text(
        value=default_user,
        label="User",
        placeholder="abcd1234",
    )
    remote_host = mo.ui.text(
        value=default_host,
        label="Host",
        placeholder="rosa.hpc.uni-oldenburg.de",
    )
    limit = mo.ui.slider(
        5,
        50,
        step=5,
        value=min(default_limit, 50),
        label="Rows",
        show_value=True,
    )
    refresh = mo.ui.run_button(label="Refresh")
    header = mo.md(f"## `{project_dir.name}`")
    advanced = mo.accordion(
        {
            "Options": mo.hstack(
                [user, remote_host, limit],
                justify="start",
                align="end",
                widths=[1, 2, 1],
            )
        }
    )
    controls_row = mo.hstack(
        [
            header,
            controller_log,
            refresh,
            advanced,
        ],
        align="end",
        justify="start",
        widths=[1, 2, 1, 2],
    )
    controls_row
    return controller_log, limit, refresh, remote_host, user


@app.cell
def _(mo):
    get_selected_details, set_selected_details = mo.state(None)
    return get_selected_details, set_selected_details


@app.cell
def _(controller_log, set_selected_details):
    controller_log.value
    set_selected_details(None)
    return


@app.cell
def _(
    apply_squeue_state,
    Path,
    controller_log,
    fetch_squeue_output,
    parse_controller_log,
    refresh,
    remote_host,
    user,
):
    refresh.value
    selected_log_path = controller_log.value or ""
    run = None
    queue_warning = ""

    if selected_log_path:
        run = parse_controller_log(Path(selected_log_path))

        queue_user = user.value.strip()
        queue_host = remote_host.value.strip()
        if queue_user and queue_host:
            squeue_output = fetch_squeue_output(user=queue_user, host=queue_host)
            if squeue_output:
                apply_squeue_state(run, squeue_output)
            else:
                queue_warning = "Could not fetch live SLURM state."
        elif queue_user and not queue_host:
            queue_warning = "Set an SSH host to enrich jobs with live SLURM state."

    return queue_warning, run


@app.cell
def _(mo, queue_warning, run):
    if run is None:
        summary_content = mo.callout(
            "No controller logs found under `log/slurm`.",
            kind="warn",
        )
    else:
        counts = run.counts()
        progress = "-"
        if run.progress_total:
            progress = f"{run.progress_done}/{run.progress_total} ({run.progress_percent}%)"

        summary = mo.hstack(
            [
                mo.stat(run.state, label="State", caption=run.started_at, bordered=True),
                mo.stat(progress, label="Progress", caption=run.status_line, bordered=True),
                mo.stat(counts["running"], label="Running", bordered=True),
                mo.stat(
                    counts["submitted"] - counts["running"],
                    label="Pending",
                    bordered=True,
                ),
                mo.stat(counts["failed"], label="Failed", bordered=True),
            ],
            gap=1,
            widths="equal",
        )

        warnings = []
        if queue_warning:
            warnings.append(mo.callout(queue_warning, kind="warn"))
        if run.errors:
            errors_md = mo.md("\n".join(f"- {message}" for message in run.errors[-5:]))
            warnings.append(mo.callout(errors_md, kind="danger"))

        summary_content = mo.vstack([summary, *warnings], gap=1)

    summary_content
    return


@app.cell
def _(job_rows, limit, mo, run, set_selected_details):
    if run is None:
        jobs_panel = mo.md("")
    else:
        rows = job_rows(run, limit=limit.value)
        if not rows:
            jobs_panel = mo.callout(
                "No jobs were parsed from the selected controller log.",
                kind="info",
            )
        else:
            def make_open_button(payload):
                def _open(_value):
                    set_selected_details(payload)
                    return payload["smk_jobid"]

                return mo.ui.button(
                    label="View",
                    value=None,
                    on_click=_open,
                    tooltip="Open logs in the viewer",
                )

            visible = []
            for row in rows:
                visible.append(
                    {
                        "open": make_open_button(row),
                        "smk_jobid": row["smk_jobid"],
                        "state": row["state"],
                        "rule": row["rule"],
                        "wildcards": row["wildcards"],
                        "elapsed": row["elapsed"],
                        "node_or_reason": row["node_or_reason"],
                        "log_file": row["log_file"],
                    }
                )

            jobs_table = mo.ui.table(
                visible,
                page_size=min(limit.value, 20),
                freeze_columns_left=["open", "rule"],
                wrapped_columns=["wildcards", "node_or_reason", "log_file"],
                header_tooltip={
                    "open": "Open logs in the viewer without leaving the dashboard",
                },
                label="Jobs",
                max_height=760,
            )
            jobs_panel = mo.vstack([mo.md("### Jobs"), jobs_table], gap=1)

    return (jobs_panel,)


@app.cell
def _(
    Path,
    get_selected_details,
    mo,
    project_dir,
    read_log_tail,
    resolve_log_paths,
    set_selected_details,
):
    selected_details = get_selected_details()

    def close_viewer(_value):
        set_selected_details(None)
        return None

    if selected_details is None:
        details_panel = mo.vstack(
            [
                mo.md("### Log Viewer"),
                mo.callout("Use `View` in the jobs table to inspect logs here.", kind="info"),
            ],
            gap=1,
        )
    else:
        log_paths = resolve_log_paths(project_dir, selected_details["log"])
        tabs = {}
        for path in log_paths:
            tabs[path.name] = mo.vstack(
                [
                    mo.plain_text(str(path)),
                    mo.ui.code_editor(
                        value=read_log_tail(Path(path)),
                        language="shell",
                        disabled=True,
                        min_height=420,
                        max_height=720,
                        show_copy_button=True,
                    ),
                ],
                gap=1,
            )

        details_panel = mo.vstack(
            [
                mo.hstack(
                    [
                        mo.md(f"### {selected_details['rule']}"),
                        mo.ui.button(
                            label="Close",
                            value=None,
                            on_click=close_viewer,
                            tooltip="Hide the log viewer",
                        ),
                    ],
                    justify="space-between",
                    align="center",
                    widths=[4, 1],
                ),
                mo.md(
                    f"""
                    `{selected_details["state"]}`  
                    `SMK {selected_details["smk_jobid"]}` `SLURM {selected_details["slurm_jobid"]}`
                    """
                ),
                mo.plain_text(selected_details["wildcards"]),
                mo.tabs(tabs) if tabs else mo.callout("No readable log files for this job.", kind="warn"),
            ],
            gap=1,
        )

    return (details_panel,)


@app.cell
def _(details_panel, jobs_panel, mo):
    layout = mo.hstack(
        [jobs_panel, details_panel],
        widths=[3, 2],
        align="start",
        gap=1,
    )
    layout
    return


if __name__ == "__main__":
    app.run()
