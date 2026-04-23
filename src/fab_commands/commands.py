"""Operational fab command implementations."""

from __future__ import annotations

import glob
import os

from . import project


def fixvideo(c, pathname, overwrite=False):
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
            fixvideo(c, pathname, overwrite)
        except Exception as exc:
            print(exc)


def zarr_to_h5(zarr_folder, delete_zarr=False):
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
            zarr_to_h5(pathname, delete_zarr=delete_zarr)
        except Exception as exc:
            print(exc)


def status(c):
    project.update_globals(c)
    c.run("sacct --format='JobID,JobName,Elapsed,MaxVMSize,MaxRSS,CPUTime,State'")


def queue(c):
    project.update_globals(c)
    c.run(f"squeue -u {project.connection_user(c)}")


def create_envs(c):
    project.update_globals(c)
    c.run(
        "umask a=rwx;"
        "module load Miniforge3;"
        "export PYTHONDONTWRITEBYTECODE=1;"
        f"source activate {project.REMOTE_CONDA_ENV};"
        f"cd {project.FOLDER};"
        "CONDA_OVERRIDE_CUDA=11.2 snakemake "
        f"--workflow-profile {project.WORKFLOW_PROFILE} "
        "--executor local "
        "--cores 1 "
        "--jobs 1 "
        "--conda-create-envs-only;"
    )


def submit(c, target=""):
    project.update_globals(c)
    controller_cmd = project.snakemake_submit_command(target=target)
    launch_cmd = project.controller_launch_command(controller_cmd)
    c.run(launch_cmd)


def cancel(c, job_id):
    project.update_globals(c)
    c.run(f"scancel {job_id}")


def unlock(c):
    project.update_globals(c)
    c.run(
        "umask a=rwx;"
        "module load Miniforge3;"
        "export PYTHONDONTWRITEBYTECODE=1;"
        f"source activate {project.REMOTE_CONDA_ENV};"
        f"cd {project.FOLDER};"
        f"snakemake --workflow-profile {project.WORKFLOW_PROFILE} --executor local --cores 1 --unlock"
    )


def annotate(c):  # noqa: ARG001
    import etho_annotator.app

    etho_annotator.app.main()
