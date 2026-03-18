"""Global fab commands compatible with Fabric 3.x."""

import glob
import logging
import os
from fabric import Connection, task

logging.basicConfig(level=logging.INFO)


NAME = os.path.basename(os.getcwd())
REMOTE = ""
SERVER = ""
FOLDER = ""

WORKFLOW_PROFILE = "../snakemake-workflows/profiles/rosa"
DEFAULT_HOSTS = [{"host": "rosa.hpc.uni-oldenburg.de"}]


def _connection_user(c):
    return getattr(c, "user", None) or os.environ.get("USER", "")


def _with_user(c, user):
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
        user = _connection_user(c)
        SERVER = f"{user}@{host}:"
        REMOTE = SERVER + FOLDER
    else:
        SERVER = ""
        REMOTE = ""


def configure_project(project_name):
    update_globals(name=project_name)


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
    c.run(f"squeue -u {_connection_user(c)}")


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


def _simple_task(func, *, name, hosts=None):
    kwargs = {"name": name}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, user=""):
        func(_with_user(c, user) if user else c)

    return wrapped


def _submit_task(*, hosts=None):
    kwargs = {"name": "submit"}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, target="", user=""):
        submit(_with_user(c, user) if user else c, target=target)

    return wrapped


def _cancel_task(*, hosts=None):
    kwargs = {"name": "cancel"}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, job_id, user=""):
        cancel(_with_user(c, user) if user else c, job_id=job_id)

    return wrapped


def _fixvideos_task(pathname_pattern):
    @task(name="fixvideos")
    def wrapped(c, pathname_pattern=pathname_pattern, overwrite=False):
        fixvideos(c, pathname_pattern=pathname_pattern, overwrite=overwrite)

    return wrapped


def register_tasks(
    namespace,
    *,
    project_name,
    include=(),
    fixvideos_pattern="dat/**/*.avi",
    annotate_hosts=None,
):
    configure_project(project_name)

    task_names = set(include)
    remote_hosts = DEFAULT_HOSTS

    remote_tasks = {
        "queue": _simple_task(queue, name="queue", hosts=remote_hosts),
        "create_envs": _simple_task(
            create_envs, name="create_envs", hosts=remote_hosts
        ),
        "status": _simple_task(status, name="status", hosts=remote_hosts),
        "unlock": _simple_task(unlock, name="unlock", hosts=remote_hosts),
        "submit": _submit_task(hosts=remote_hosts),
        "cancel": _cancel_task(hosts=remote_hosts),
    }
    local_tasks = {
        "annotate": _simple_task(annotate, name="annotate", hosts=annotate_hosts),
        "fixvideos": _fixvideos_task(fixvideos_pattern),
    }

    for task_name in task_names:
        if task_name in remote_tasks:
            namespace[task_name] = remote_tasks[task_name]
        elif task_name in local_tasks:
            namespace[task_name] = local_tasks[task_name]
        else:
            raise ValueError(f"Unknown fab task: {task_name}")
