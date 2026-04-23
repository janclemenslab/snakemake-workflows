"""Fabric task wrappers for terminal Snakemake monitoring."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

from fabric import task

from .monitor import apply_squeue_state, latest_controller_log, parse_controller_log, render_monitor


ProjectDirGetter = Callable[[], Path]
ConnectionUserGetter = Callable[[Any], str]
ConnectionWrapper = Callable[[Any, str], Any]


def run_monitor(
    c: Any,
    *,
    project_dir: ProjectDirGetter,
    connection_user: ConnectionUserGetter,
    limit: int | str = 20,
) -> None:
    workflow_dir = project_dir()
    controller_log = latest_controller_log(workflow_dir)
    if controller_log is None:
        print("No controller logs found under log/slurm.")
        return

    run = parse_controller_log(controller_log)

    if run.needs_live_slurm_state():
        try:
            result = c.run(
                f"squeue -h -u {connection_user(c)} -o '%i|%T|%M|%N|%R'",
                hide=True,
                warn=True,
            )
            if result.ok:
                apply_squeue_state(run, result.stdout)
        except Exception:
            pass

    render_monitor(run, workflow_name=workflow_dir.name, limit=int(limit))


def build_monitor_task(
    *,
    project_dir: ProjectDirGetter,
    connection_user: ConnectionUserGetter,
    with_user: ConnectionWrapper,
    hosts=None,
):
    kwargs = {"name": "monitor"}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, limit=20, user=""):
        run_monitor(
            with_user(c, user) if user else c,
            project_dir=project_dir,
            connection_user=connection_user,
            limit=limit,
        )

    return wrapped
