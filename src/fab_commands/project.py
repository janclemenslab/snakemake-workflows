"""Shared project state and command-building helpers for fab commands."""

from __future__ import annotations

import os
from pathlib import Path
import shlex

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


def snakemake_target_arg(target):
    target = (target or "").strip()
    return f" {shlex.quote(target)}" if target else ""


def snakemake_submit_command(target=""):
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
        f"{snakemake_target_arg(target)}"
    )
