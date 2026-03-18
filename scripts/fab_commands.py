"""Global fab commands compatible with Fabric 3.x."""

import glob
import logging
import os
from pathlib import Path

logging.basicConfig(level=logging.INFO)


# define globals
NAME = os.path.basename(os.getcwd())
REMOTE = ""
SERVER = ""
FOLDER = ""

WORKFLOW_PROFILE = "../snakemake-workflows/profiles/rosa"


def update_globals(c=None, name=None):
    """Register derived globals from project name and optional connection."""
    global SERVER, REMOTE, FOLDER, NAME
    if name:
        NAME = name

    FOLDER = f"/fs/s6k/groups/agauneu/#Data/{NAME}/"

    if c is not None and getattr(c, "host", None):
        host = getattr(c, "original_host", None) or c.host
        user = getattr(c, "user", None) or os.environ.get("USER", "")
        SERVER = f"{user}@{host}:"
        REMOTE = SERVER + FOLDER
    else:
        SERVER = ""
        REMOTE = ""


def _snakemake_target_arg(target):
    target = (target or "").strip()
    return f" {target}" if target else ""


def _fixvideo(c, pathname, overwrite=False):
    """Convert AVI to MP4 using stream copy."""
    trunk = os.path.splitext(pathname)[0]
    if not overwrite and not os.path.exists(f"{trunk}.mp4"):
        c.local(f'ffmpeg -i "{pathname}" -c copy {trunk}.mp4', warn=True)


def fixvideos(c, pathname_pattern="dat/**/*.avi", overwrite=False):
    """Batch-convert AVI files under dat/ to MP4."""
    pathnames = glob.glob(pathname_pattern, recursive=True)
    print(pathnames)
    for pathname in pathnames:
        try:
            _fixvideo(c, pathname, overwrite)
        except Exception as exc:
            print(exc)


def status(c):
    update_globals(c)
    c.run("sacct --format='JobID,JobName,Elapsed,MaxVMSize,MaxRSS,CPUTime,State'")


def queue(c):
    update_globals(c)
    c.run(f"squeue -u {c.user}")


def create_envs(c):
    update_globals(c)
    c.run(
        "umask a=rwx;"
        "module load Miniforge3;"
        "export PYTHONDONTWRITEBYTECODE=1;"
        "source activate /fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb;"
        f"cd {FOLDER};"
        "CONDA_OVERRIDE_CUDA=11.2 snakemake "
        f"--workflow-profile {WORKFLOW_PROFILE} "
        "--executor local "
        "--cores 1 "
        "--jobs 1 "
        "--conda-create-envs-only;"
    )


def submit(c, target=""):
    update_globals(c)
    c.run(
        "umask a=rwx;"
        "module load Miniforge3;"
        "export PYTHONDONTWRITEBYTECODE=1;"
        "source activate /fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb;"
        f"cd {FOLDER};"
        "CONDA_OVERRIDE_CUDA=11.2 snakemake "
        f"--workflow-profile {WORKFLOW_PROFILE} "
        "--conda-create-envs-only;"
        "CONDA_OVERRIDE_CUDA=11.2 snakemake "
        f"--workflow-profile {WORKFLOW_PROFILE}"
        f"{_snakemake_target_arg(target)}"
    )


def cancel(c, job_id):
    update_globals(c)
    c.run(f"scancel {job_id}")


def unlock(c):
    update_globals(c)
    c.run(
        "umask a=rwx;"
        "module load Miniforge3;"
        "export PYTHONDONTWRITEBYTECODE=1;"
        "source activate /fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb;"
        f"cd {FOLDER};"
        f"snakemake --workflow-profile {WORKFLOW_PROFILE} --executor local --cores 1 --unlock"
    )


def annotate(c):  # noqa: ARG001
    import etho_annotator.app

    etho_annotator.app.main()
