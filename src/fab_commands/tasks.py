"""Fabric task registration for shared fab commands."""

from __future__ import annotations

from fabric import task

from . import commands, project
from .dashboard_task import build_dashboard_task
from .monitor_task import build_monitor_task


def simple_task(func, *, name, hosts=None):
    kwargs = {"name": name}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, user=""):
        func(project.with_user(c, user) if user else c)

    return wrapped


def submit_task(*, hosts=None):
    kwargs = {"name": "submit"}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, target="", user=""):
        commands.submit(project.with_user(c, user) if user else c, target=target)

    return wrapped


def cancel_task(*, hosts=None):
    kwargs = {"name": "cancel"}
    if hosts is not None:
        kwargs["hosts"] = hosts

    @task(**kwargs)
    def wrapped(c, job_id, user=""):
        commands.cancel(project.with_user(c, user) if user else c, job_id=job_id)

    return wrapped


def fixvideos_task(pathname_pattern):
    @task(name="fixvideos")
    def wrapped(c, pathname_pattern=pathname_pattern, overwrite=False):
        commands.fixvideos(c, pathname_pattern=pathname_pattern, overwrite=overwrite)

    return wrapped


def fixdaq_task(pathname_pattern):
    @task(name="fixdaq")
    def wrapped(c, pathname_pattern=pathname_pattern, delete_zarr=True):
        commands.fixdaq(c, pathname_pattern=pathname_pattern, delete_zarr=delete_zarr)

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
    project.configure_project(project_name, project_dir=project_dir)

    task_names = set(include)
    remote_hosts = project.DEFAULT_HOSTS

    remote_tasks = {
        "queue": simple_task(commands.queue, name="queue", hosts=remote_hosts),
        "monitor": build_monitor_task(
            project_dir=project.local_project_dir,
            connection_user=project.connection_user,
            with_user=project.with_user,
            hosts=remote_hosts,
        ),
        "create_envs": simple_task(
            commands.create_envs, name="create_envs", hosts=remote_hosts
        ),
        "status": simple_task(commands.status, name="status", hosts=remote_hosts),
        "unlock": simple_task(commands.unlock, name="unlock", hosts=remote_hosts),
        "submit": submit_task(hosts=remote_hosts),
        "cancel": cancel_task(hosts=remote_hosts),
    }
    local_tasks = {
        "annotate": simple_task(
            commands.annotate, name="annotate", hosts=annotate_hosts
        ),
        "dashboard": build_dashboard_task(
            project_dir=project.local_project_dir,
            default_hosts=project.DEFAULT_HOSTS,
        ),
        "fixvideos": fixvideos_task(fixvideos_pattern),
        "fixdaq": fixdaq_task("dat/**/*daq.zarr"),
    }

    for task_name in task_names:
        if task_name in remote_tasks:
            namespace[task_name] = remote_tasks[task_name]
        elif task_name in local_tasks:
            namespace[task_name] = local_tasks[task_name]
        else:
            raise ValueError(f"Unknown fab task: {task_name}")
