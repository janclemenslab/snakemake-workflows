"""Shared project state and command-building helpers for fab commands."""

from __future__ import annotations

import os
from pathlib import Path
import shlex
from typing import Iterable

from fabric import Connection


NAME = os.path.basename(os.getcwd())
REMOTE = ""
SERVER = ""
FOLDER = ""
LOCAL_PROJECT_DIR = Path(os.getcwd()).resolve()

WORKFLOW_PROFILE = "../snakemake-workflows/profiles/rosa"
REMOTE_CONDA_ENV = "/fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb"
CONTROLLER_LOG_DIR = "log/slurm"
DEFAULT_HOSTS = [{"host": "rosa.hpc.uni-oldenburg.de"}]


def connection_user(c):
    return getattr(c, "user", None) or os.environ.get("USER", "")


def with_user(c, user):
    host = getattr(c, "original_host", None) or getattr(c, "host", None)
    if not host:
        return c
    return Connection(host=host, user=user, config=c.config)


def update_globals(c=None, name=None):
    """Register derived globals from project name and optional connection."""
    global SERVER, REMOTE, FOLDER, NAME
    if name:
        NAME = name

    FOLDER = f"/fs/s6k/groups/agauneu/#Data/{NAME}/"

    if c is not None and getattr(c, "host", None):
        host = getattr(c, "original_host", None) or c.host
        user = connection_user(c)
        SERVER = f"{user}@{host}:"
        REMOTE = SERVER + FOLDER
    else:
        SERVER = ""
        REMOTE = ""


def configure_project(project_name, project_dir=None):
    global LOCAL_PROJECT_DIR
    update_globals(name=project_name)
    if project_dir is not None:
        LOCAL_PROJECT_DIR = Path(project_dir).resolve()


def local_project_dir() -> Path:
    return LOCAL_PROJECT_DIR


def _normalize_shell_args(args: str | Iterable[str] | None) -> list[str]:
    if args is None:
        return []
    if isinstance(args, str):
        return [arg for arg in shlex.split(args) if arg.strip()]
    return [str(arg).strip() for arg in args if str(arg).strip()]


def snakemake_target_args(
    target: str | Iterable[str] | None = "",
    *,
    targets: str | Iterable[str] | None = None,
) -> str:
    effective_targets = _normalize_shell_args(targets if targets is not None else target)
    return "".join(f" {shlex.quote(value)}" for value in effective_targets)


def snakemake_submit_command(
    target: str | Iterable[str] | None = "",
    *,
    targets: str | Iterable[str] | None = None,
    force: bool = False,
) -> str:
    target_args = snakemake_target_args(target, targets=targets)
    force_args = " --force" if force and target_args else ""
    return (
        "umask a=rwx;"
        "module load Miniforge3;"
        "export PYTHONDONTWRITEBYTECODE=1;"
        f"source activate {shlex.quote(REMOTE_CONDA_ENV)};"
        f"cd {shlex.quote(FOLDER)};"
        "CONDA_OVERRIDE_CUDA=11.2 snakemake "
        f"--workflow-profile {shlex.quote(WORKFLOW_PROFILE)} "
        "--conda-create-envs-only;"
        "CONDA_OVERRIDE_CUDA=11.2 snakemake "
        f"--workflow-profile {shlex.quote(WORKFLOW_PROFILE)}"
        f"{force_args}{target_args}"
    )


def controller_launch_command(controller_cmd: str) -> str:
    return "bash -lc " + shlex.quote(
        f"cd {shlex.quote(FOLDER)};"
        f"mkdir -p {shlex.quote(CONTROLLER_LOG_DIR)};"
        'controller_log="log/slurm/controller-$(date +%Y%m%dT%H%M%S).log";'
        f"nohup bash -lc {shlex.quote(controller_cmd)} "
        '> "$controller_log" 2>&1 < /dev/null &'
        "controller_pid=$!;"
        'printf "Started Snakemake controller PID %s\\nLog: %s\\n" '
        '"$controller_pid" "$controller_log"'
    )
