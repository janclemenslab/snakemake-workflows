"""Fabric task wrappers for the Snakemake dashboard."""

from __future__ import annotations

from collections.abc import Callable, Sequence
import logging
from pathlib import Path
import shlex
from typing import Any

from fabric import task


ProjectDirGetter = Callable[[], Path]


def _browser_url(bind_host: str, port: int) -> str:
    host = (bind_host or "").strip()
    if not host or host in {"0.0.0.0", "::"}:
        host = "127.0.0.1"
    return f"http://{host}:{port}/"


def launch_dashboard(
    c: Any,
    *,
    project_dir: ProjectDirGetter,
    default_hosts: Sequence[dict[str, str]] | None = None,
    user: str = "",
    port: int = 2718,
    bind_host: str = "127.0.0.1",
    remote_host: str = "",
    limit: int = 100,
    headless: bool = False,
    open_browser: bool = True,
    verbose: bool = False,
) -> None:
    dashboard_app = Path(__file__).with_name("monitor_dashboard.py")
    workflow_dir = project_dir()
    effective_remote_host = remote_host or (
        default_hosts[0]["host"] if default_hosts else ""
    )
    effective_open_browser = bool(open_browser) and not bool(headless)
    command = [
        "env",
        f"FAB_MONITOR_PROJECT_DIR={workflow_dir}",
        f"FAB_MONITOR_USER={user}",
        f"FAB_MONITOR_REMOTE_HOST={effective_remote_host}",
        f"FAB_MONITOR_LIMIT={limit}",
        f"FAB_MONITOR_BIND_HOST={bind_host}",
        f"FAB_MONITOR_PORT={port}",
        f"FAB_MONITOR_OPEN_BROWSER={1 if effective_open_browser else 0}",
        f"FAB_MONITOR_BROWSER_URL={_browser_url(bind_host, int(port))}",
        f"FAB_MONITOR_VERBOSE={1 if verbose else 0}",
        "python",
        str(dashboard_app),
    ]
    if headless:
        logging.info(
            "Suppressing browser launch because headless=%s",
            headless,
        )

    c.run(" ".join(shlex.quote(part) for part in command))


def build_dashboard_task(
    *,
    project_dir: ProjectDirGetter,
    default_hosts: Sequence[dict[str, str]] | None = None,
):
    @task(name="dashboard")
    def wrapped(
        c,
        user="",
        port=2718,
        bind_host="127.0.0.1",
        remote_host="",
        limit=100,
        headless=False,
        open_browser=True,
        verbose=False,
    ):
        launch_dashboard(
            c,
            project_dir=project_dir,
            default_hosts=default_hosts,
            user=user,
            port=port,
            bind_host=bind_host,
            remote_host=remote_host,
            limit=limit,
            headless=headless,
            open_browser=open_browser,
            verbose=verbose,
        )

    return wrapped
