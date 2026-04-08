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
    import pandas as pd

    from fab_commands.monitor import (
        apply_squeue_state,
        controller_logs,
        fetch_squeue_output,
        job_rows,
        parse_controller_log,
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
        pd,
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
        label="Controller log",
        full_width=True,
    )
    user = mo.ui.text(
        value=default_user,
        label="SLURM user",
        placeholder="abcd1234",
        full_width=True,
    )
    remote_host = mo.ui.text(
        value=default_host,
        label="SSH host",
        placeholder="rosa.hpc.uni-oldenburg.de",
        full_width=True,
    )
    limit = mo.ui.slider(
        5,
        100,
        step=5,
        value=min(default_limit, 100),
        label="Visible jobs",
        show_value=True,
    )
    refresh = mo.ui.run_button(label="Refresh now")
    header = mo.md(
        f"""
        # Snakemake dashboard

        **Workflow:** `{project_dir.name}`  
        **Project dir:** `{project_dir}`
        """
    )
    controls = mo.vstack(
        [
            header,
            mo.hstack([controller_log, user, remote_host], align="end", widths=[3, 1, 2]),
            mo.hstack([limit, refresh], justify="start", align="center", widths=[3, 1]),
        ],
        gap=1,
    )
    controls
    return controller_log, limit, refresh, remote_host, user


@app.cell
def _(
    apply_squeue_state,
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

        summary = mo.md(
            f"""
            ## Summary

            **State:** `{run.state}`  
            **Started:** `{run.started_at}`  
            **Progress:** `{progress}`  
            **Last event:** `{run.status_line}`  
            **Counts:** `running={counts["running"]}` `pending={counts["submitted"] - counts["running"]}` `finished={counts["finished"]}` `failed={counts["failed"]}` `local={counts["local"]}`
            """
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
def _(job_rows, limit, mo, pd, run):
    if run is None:
        jobs_table = None
        jobs_content = mo.md("")
    else:
        rows = job_rows(run, limit=limit.value)
        frame = pd.DataFrame(rows)
        if frame.empty:
            jobs_table = None
            jobs_content = mo.callout(
                "No jobs were parsed from the selected controller log.",
                kind="info",
            )
        else:
            jobs_table = mo.ui.table(
                frame,
                page_size=min(limit.value, 20),
                selection="single",
                wrapped_columns=["wildcards", "node_or_reason", "log"],
                label="Jobs",
            )
            jobs_content = jobs_table

    jobs_content
    return (jobs_table,)


@app.cell
def _(jobs_table, mo):
    selected_row = None
    if jobs_table is not None and jobs_table.value:
        selected_row = jobs_table.value
        if isinstance(selected_row, list):
            selected_row = selected_row[0] if selected_row else None

    if selected_row is None:
        details = mo.md("")
    else:
        details = mo.md(
            f"""
            ## Selected job

            **Rule:** `{selected_row["rule"]}`  
            **Snakemake job:** `{selected_row["smk_jobid"]}`  
            **SLURM job:** `{selected_row["slurm_jobid"]}`  
            **State:** `{selected_row["state"]}`  
            **Wildcards:** `{selected_row["wildcards"]}`  
            **Input:** `{selected_row["input"]}`  
            **Output:** `{selected_row["output"]}`  
            **Log:** `{selected_row["log"]}`
            """
        )

    details
    return


if __name__ == "__main__":
    app.run()
