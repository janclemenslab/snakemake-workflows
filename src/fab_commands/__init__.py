"""Global fab commands compatible with Fabric 3.x."""

import glob
import logging
import os
from pathlib import Path
import shlex
from fabric import Connection, task

from .monitor import apply_squeue_state, latest_controller_log, parse_controller_log, render_monitor

logging.basicConfig(level=logging.INFO)


NAME = os.path.basename(os.getcwd())
REMOTE = ""
SERVER = ""
FOLDER = ""
LOCAL_PROJECT_DIR = Path(os.getcwd()).resolve()

WORKFLOW_PROFILE = "../snakemake-workflows/profiles/rosa"
REMOTE_CONDA_ENV = "/fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb"
CONTROLLER_LOG_DIR = "log/slurm"
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


def configure_project(project_name, project_dir=None):
    global LOCAL_PROJECT_DIR
    update_globals(name=project_name)
    if project_dir is not None:
        LOCAL_PROJECT_DIR = Path(project_dir).resolve()


def _local_project_dir() -> Path:
    return LOCAL_PROJECT_DIR


def _snakemake_target_arg(target):
    target = (target or "").strip()
    return f" {shlex.quote(target)}" if target else ""


def _snakemake_submit_command(target=""):
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
        f"{_snakemake_target_arg(target)}"
    )


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


def _zarr_to_h5(zarr_folder, delete_zarr=False):
    """Convert a Zarr folder to HDF5, copying array attributes."""
    import shutil

    import h5py
    import numpy as np

    try:
        import zarr
    except ImportError as exc:
        raise RuntimeError(
            "fixdaq requires the optional 'zarr' package in the active environment"
        ) from exc

    def copy_attrs(source_attrs, target_attrs):
        for key, value in source_attrs.items():
            if isinstance(value, list):
                if all(isinstance(item, str) for item in value):
                    target_attrs[key] = np.array(
                        value, dtype=h5py.string_dtype(encoding="utf-8")
                    )
                else:
                    target_attrs[key] = np.array(value)
            elif value is None:
                target_attrs[key] = "None"
            else:
                target_attrs[key] = value

    def copy_group(zgroup, h5group):
        copy_attrs(zgroup.attrs, h5group.attrs)
        for name, item in zgroup.items():
            if isinstance(item, zarr.Array):
                print(f"Copying array: {name}")
                dset = h5group.create_dataset(name, data=item[:])
                copy_attrs(item.attrs, dset.attrs)
            elif isinstance(item, zarr.Group):
                new_group = h5group.create_group(name)
                copy_group(item, new_group)

    h5_path = os.path.splitext(zarr_folder)[0] + ".h5"
    if os.path.exists(zarr_folder) and not os.path.exists(h5_path):
        root = zarr.open(zarr_folder, mode="r")
        with h5py.File(h5_path, "w") as h5_file:
            copy_group(root, h5_file)
        if delete_zarr:
            try:
                with h5py.File(h5_path, "r") as h5_file:
                    len(h5_file.keys())
                print("Deleting zarr folder")
                shutil.rmtree(zarr_folder)
            except Exception:
                pass
        print(f"Conversion complete: {h5_path}")


def fixdaq(c, pathname_pattern="dat/**/*daq.zarr", delete_zarr=True):  # noqa: ARG001
    """Batch-convert DAQ Zarr folders under dat/ to HDF5."""
    pathnames = glob.glob(pathname_pattern, recursive=True)
    print(pathnames)
    for pathname in pathnames:
        try:
            _zarr_to_h5(pathname, delete_zarr=delete_zarr)
        except Exception as exc:
            print(exc)


def status(c):
    update_globals(c)
    c.run("sacct --format='JobID,JobName,Elapsed,MaxVMSize,MaxRSS,CPUTime,State'")


def monitor(c, limit=20):
    project_dir = _local_project_dir()
    controller_log = latest_controller_log(project_dir)
    if controller_log is None:
        print("No controller logs found under log/slurm.")
        return

    run = parse_controller_log(controller_log)

    try:
        result = c.run(
            f"squeue -h -u {_connection_user(c)} -o '%i|%T|%M|%N|%R'",
            hide=True,
            warn=True,
        )
        if result.ok:
            apply_squeue_state(run, result.stdout)
    except Exception:
        pass

    render_monitor(run, workflow_name=project_dir.name, limit=int(limit))


def dashboard(
    c,
    user="",
    port=2718,
    bind_host="127.0.0.1",
    remote_host="",
    limit=25,
    headless=False,
):
    dashboard_app = Path(__file__).with_name("monitor_dashboard.py")
    project_dir = _local_project_dir()
    effective_remote_host = remote_host or (
        DEFAULT_HOSTS[0]["host"] if DEFAULT_HOSTS else ""
    )
    command = [
        "env",
        f"FAB_MONITOR_PROJECT_DIR={project_dir}",
        f"FAB_MONITOR_USER={user}",
        f"FAB_MONITOR_REMOTE_HOST={effective_remote_host}",
        f"FAB_MONITOR_LIMIT={limit}",
        "marimo",
        "run",
        str(dashboard_app),
        "--host",
        str(bind_host),
        "--port",
        str(port),
    ]
    if headless:
        command.append("--headless")
    command.extend(
        [
            "--",
            "--project-dir",
            str(project_dir),
            "--limit",
            str(limit),
        ]
    )
    if user:
        command.extend(["--user", user])
    if effective_remote_host:
        command.extend(["--remote-host", effective_remote_host])

    c.run(" ".join(shlex.quote(part) for part in command))


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
    controller_cmd = _snakemake_submit_command(target=target)
    launch_cmd = "bash -lc " + shlex.quote(
        f"cd {shlex.quote(FOLDER)};"
        f"mkdir -p {shlex.quote(CONTROLLER_LOG_DIR)};"
        'controller_log="log/slurm/controller-$(date +%Y%m%dT%H%M%S).log";'
        f"nohup bash -lc {shlex.quote(controller_cmd)} "
        '> "$controller_log" 2>&1 < /dev/null &'
        "controller_pid=$!;"
        'printf "Started Snakemake controller PID %s\\nLog: %s\\n" '
        '"$controller_pid" "$controller_log"'
    )
    c.run(launch_cmd)


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


def _monitor_task(*, hosts=None):
    kwargs = {"name": "monitor"}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, limit=20, user=""):
        monitor(_with_user(c, user) if user else c, limit=limit)

    return wrapped


def _cancel_task(*, hosts=None):
    kwargs = {"name": "cancel"}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, job_id, user=""):
        cancel(_with_user(c, user) if user else c, job_id=job_id)

    return wrapped


def _dashboard_task():
    @task(name="dashboard")
    def wrapped(
        c,
        user="",
        port=2718,
        bind_host="127.0.0.1",
        remote_host="",
        limit=25,
        headless=False,
    ):
        dashboard(
            c,
            user=user,
            port=port,
            bind_host=bind_host,
            remote_host=remote_host,
            limit=limit,
            headless=headless,
        )

    return wrapped


def _fixvideos_task(pathname_pattern):
    @task(name="fixvideos")
    def wrapped(c, pathname_pattern=pathname_pattern, overwrite=False):
        fixvideos(c, pathname_pattern=pathname_pattern, overwrite=overwrite)

    return wrapped


def _fixdaq_task(pathname_pattern):
    @task(name="fixdaq")
    def wrapped(c, pathname_pattern=pathname_pattern, delete_zarr=True):
        fixdaq(c, pathname_pattern=pathname_pattern, delete_zarr=delete_zarr)

    return wrapped


def register_tasks(
    namespace,
    *,
    project_name,
    project_dir=None,
    include=(),
    fixvideos_pattern="dat/**/*.avi",
    annotate_hosts=None,
):
    configure_project(project_name, project_dir=project_dir)

    task_names = set(include)
    remote_hosts = DEFAULT_HOSTS

    remote_tasks = {
        "queue": _simple_task(queue, name="queue", hosts=remote_hosts),
        "monitor": _monitor_task(hosts=remote_hosts),
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
        "dashboard": _dashboard_task(),
        "fixvideos": _fixvideos_task(fixvideos_pattern),
        "fixdaq": _fixdaq_task("dat/**/*daq.zarr"),
    }

    for task_name in task_names:
        if task_name in remote_tasks:
            namespace[task_name] = remote_tasks[task_name]
        elif task_name in local_tasks:
            namespace[task_name] = local_tasks[task_name]
        else:
            raise ValueError(f"Unknown fab task: {task_name}")
