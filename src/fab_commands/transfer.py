"""Transfer helpers for moving local data into the shared SMB tree."""

from __future__ import annotations

import argparse
import glob
import logging
import os
import shlex
import shutil
from datetime import datetime
from pathlib import Path, PurePosixPath

from fabric import Connection


logger = logging.getLogger(__name__)

HOME_TRANSFER_DIR = Path.home() / "data.transfer"
FIXDAQ_PATTERN = "dat/**/*daq.zarr"
DEFAULT_RIG_NAME = "chainingmic"
DEFAULT_PERMISSION_HOST = "rosa.hpc.uni-oldenburg.de"
DEFAULT_REMOTE_DATA_ROOT = PurePosixPath("/fs/s6k/groups/agauneu/#Data")
DEFAULT_SMB_DATA_ROOT = Path(r"\\dss.hpc.uni-oldenburg.de\groups\agauneu\#Data")
EXCLUDED_DIRS = {"logs"}

_LOGFILE_PATH: Path | None = None


def ensure_logging_configured():
    global _LOGFILE_PATH
    if _LOGFILE_PATH is not None:
        return _LOGFILE_PATH

    logs_dir = HOME_TRANSFER_DIR / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    logfile_path = logs_dir / f"{datetime.now().strftime('%Y%m%d-%H%M%S')}.log"

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(logfile_path, encoding="utf-8")
    file_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    _LOGFILE_PATH = logfile_path
    logger.info("Logging to: %s", logfile_path)
    return logfile_path


def _copy_attrs(source_attrs, target_attrs):
    import h5py
    import numpy as np

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


def _zarr_to_h5(zarr_folder, delete_zarr=False):
    """Convert a Zarr file to HDF5, keeping attributes as dataset attributes."""
    import zarr

    ensure_logging_configured()
    logger.info("Converting zarr to h5: %s", zarr_folder)

    def copy_group(zgroup, h5group):
        _copy_attrs(zgroup.attrs, h5group.attrs)
        for name, item in zgroup.items():
            if isinstance(item, zarr.Array):
                logger.info("Copying array: %s", name)
                dset = h5group.create_dataset(name, data=item[:])
                _copy_attrs(item.attrs, dset.attrs)
            elif isinstance(item, zarr.Group):
                new_group = h5group.create_group(name)
                copy_group(item, new_group)

    h5_path = os.path.splitext(zarr_folder)[0] + ".h5"
    if os.path.exists(zarr_folder) and not os.path.exists(h5_path):
        import h5py

        root = zarr.open(zarr_folder, mode="r")
        with h5py.File(h5_path, "w") as h5_file:
            copy_group(root, h5_file)
        if delete_zarr:
            try:
                with h5py.File(h5_path, "r") as h5_file:
                    len(h5_file.keys())
                logger.info("Deleting zarr folder")
                shutil.rmtree(zarr_folder)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Could not delete zarr folder %s: %s", zarr_folder, exc)
        logger.info("Conversion complete: %s", h5_path)
    else:
        logger.info("Skipping conversion (missing source or h5 exists): %s", zarr_folder)


def transfer_fixdaq(
    source_root=HOME_TRANSFER_DIR,
    pathname_pattern=FIXDAQ_PATTERN,
    delete_zarr=True,
):
    ensure_logging_configured()
    source_root = Path(source_root)
    logger.info("Starting transfer fixdaq in %s, delete_zarr=%s", source_root, delete_zarr)

    pattern = str(source_root / pathname_pattern)
    pathnames = glob.glob(pattern, recursive=True)

    logger.info("Found %s matching path(s)", len(pathnames))
    for pathname in pathnames:
        try:
            _zarr_to_h5(pathname, delete_zarr=delete_zarr)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Failed converting %s: %s", pathname, exc)
    logger.info("Finished transfer fixdaq")


def _rig_permission_targets(rig_name, remote_data_root):
    remote_root = PurePosixPath(str(remote_data_root)) / rig_name
    return [remote_root / "dat", remote_root / "res"]


def _repair_remote_permissions(
    rig_name,
    *,
    permission_user="",
    permission_host=DEFAULT_PERMISSION_HOST,
    remote_data_root=DEFAULT_REMOTE_DATA_ROOT,
):
    ensure_logging_configured()
    connection_kwargs = {"host": permission_host}
    if permission_user:
        connection_kwargs["user"] = permission_user
    connection = Connection(**connection_kwargs)

    quoted_targets = " ".join(
        shlex.quote(target.as_posix())
        for target in _rig_permission_targets(rig_name, remote_data_root)
    )
    command = (
        "set -e;"
        f"for target in {quoted_targets}; do "
        'if [ -d "$target" ]; then '
        'chmod g+rwx "$target" >/dev/null 2>&1 || true;'
        'chmod a+rx "$target" >/dev/null 2>&1 || true;'
        "fi; "
        "done"
    )

    try:
        logger.info(
            "Repairing permissions for %s via %s",
            rig_name,
            permission_host,
        )
        result = connection.run(command, warn=True, hide=True)
        if result.stdout.strip():
            logger.info(result.stdout.strip())
        if result.stderr.strip():
            logger.warning(result.stderr.strip())
        if result.failed:
            logger.error("Permission repair failed for %s", rig_name)
            return False
        logger.info("Permissions repaired for %s", rig_name)
        return True
    except Exception as exc:  # noqa: BLE001
        logger.exception("Failed to repair permissions for %s: %s", rig_name, exc)
        return False
    finally:
        connection.close()


def transfer_fixperms(
    *,
    rig_name=DEFAULT_RIG_NAME,
    remote_data_root=DEFAULT_REMOTE_DATA_ROOT,
    permission_user="",
    permission_host=DEFAULT_PERMISSION_HOST,
):
    return _repair_remote_permissions(
        rig_name,
        permission_user=permission_user,
        permission_host=permission_host,
        remote_data_root=remote_data_root,
    )


def transfer_push(
    *,
    rig_name=DEFAULT_RIG_NAME,
    source_root=HOME_TRANSFER_DIR,
    smb_data_root=DEFAULT_SMB_DATA_ROOT,
    remote_data_root=DEFAULT_REMOTE_DATA_ROOT,
    permission_user="",
    permission_host=DEFAULT_PERMISSION_HOST,
    repair_permissions=True,
):
    ensure_logging_configured()
    logger.info("Starting transfer push")
    source_root = Path(source_root)
    destination_root = Path(smb_data_root) / rig_name / "dat"

    if not source_root.exists() or not source_root.is_dir():
        logger.error("Source folder not found: %s", source_root)
        return

    if not destination_root.exists() or not destination_root.is_dir():
        logger.error("Destination folder not accessible: %s", destination_root)
        return

    dirs_to_move = sorted(path for path in source_root.iterdir() if path.is_dir())
    if not dirs_to_move:
        logger.info("No directories to move in: %s", source_root)
        return

    for source_dir in dirs_to_move:
        if source_dir.name.lower() in EXCLUDED_DIRS:
            logger.info("Skipping %s: excluded", source_dir.name)
            continue

        destination_dir = destination_root / source_dir.name

        try:
            if destination_dir.exists():
                logger.info("Skipping %s: destination exists", source_dir.name)
            else:
                shutil.move(str(source_dir), str(destination_root))
                logger.info("Moved: %s -> %s", source_dir, destination_root)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Failed to move %s: %s", source_dir, exc)

    if repair_permissions:
        _repair_remote_permissions(
            rig_name,
            permission_user=permission_user,
            permission_host=permission_host,
            remote_data_root=remote_data_root,
        )

    logger.info("Finished transfer push")


def transfer_all(
    *,
    rig_name=DEFAULT_RIG_NAME,
    source_root=HOME_TRANSFER_DIR,
    smb_data_root=DEFAULT_SMB_DATA_ROOT,
    remote_data_root=DEFAULT_REMOTE_DATA_ROOT,
    permission_user="",
    permission_host=DEFAULT_PERMISSION_HOST,
    repair_permissions=True,
    delete_zarr=True,
):
    transfer_fixdaq(source_root=source_root, delete_zarr=delete_zarr)
    transfer_push(
        rig_name=rig_name,
        source_root=source_root,
        smb_data_root=smb_data_root,
        remote_data_root=remote_data_root,
        permission_user=permission_user,
        permission_host=permission_host,
        repair_permissions=repair_permissions,
    )


def _resolve_common_kwargs(args):
    return {
        "rig_name": args.rig_name,
        "source_root": args.source_root,
        "smb_data_root": args.smb_data_root,
        "remote_data_root": args.remote_data_root,
        "permission_user": args.permission_user,
        "permission_host": args.permission_host,
    }


def build_parser():
    parser = argparse.ArgumentParser(description="Transfer helper script")
    subparsers = parser.add_subparsers(dest="command", required=True)

    fixdaq_parser = subparsers.add_parser("fixdaq", help="Convert DAQ zarr files to h5")
    fixdaq_parser.add_argument(
        "--source-root",
        default=str(HOME_TRANSFER_DIR),
        help=f"Source transfer root (default: {HOME_TRANSFER_DIR})",
    )
    fixdaq_parser.add_argument(
        "--delete-zarr",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Delete zarr folder after successful conversion",
    )

    for command_name, help_text in {
        "push": "Move data.transfer subfolders to destination and repair permissions",
        "all": "Run fixdaq first, then push and repair permissions",
    }.items():
        command_parser = subparsers.add_parser(command_name, help=help_text)
        command_parser.add_argument(
            "--rig-name",
            default=DEFAULT_RIG_NAME,
            help=f"Destination rig folder under shared #Data (default: {DEFAULT_RIG_NAME})",
        )
        command_parser.add_argument(
            "--source-root",
            default=str(HOME_TRANSFER_DIR),
            help=f"Source transfer root (default: {HOME_TRANSFER_DIR})",
        )
        command_parser.add_argument(
            "--smb-data-root",
            default=str(DEFAULT_SMB_DATA_ROOT),
            help=f"SMB root for moved data (default: {DEFAULT_SMB_DATA_ROOT})",
        )
        command_parser.add_argument(
            "--remote-data-root",
            default=str(DEFAULT_REMOTE_DATA_ROOT),
            help=f"Server-side root used for the rig folder chmod step (default: {DEFAULT_REMOTE_DATA_ROOT})",
        )
        command_parser.add_argument(
            "--permission-user",
            default="",
            help="Optional SSH username for the permission repair host",
        )
        command_parser.add_argument(
            "--permission-host",
            default=DEFAULT_PERMISSION_HOST,
            help=f"SSH host used for permission repair (default: {DEFAULT_PERMISSION_HOST})",
        )
        command_parser.add_argument(
            "--repair-permissions",
            action=argparse.BooleanOptionalAction,
            default=True,
            help="Run remote permission repair after each move",
        )
        if command_name == "all":
            command_parser.add_argument(
                "--delete-zarr",
                action=argparse.BooleanOptionalAction,
                default=True,
                help="Delete zarr folder after successful conversion",
            )

    fixperms_parser = subparsers.add_parser(
        "fixperms",
        help="Repair permissions on the rig's dat and res folders",
    )
    fixperms_parser.add_argument(
        "--rig-name",
        default=DEFAULT_RIG_NAME,
        help=f"Destination rig folder under shared #Data (default: {DEFAULT_RIG_NAME})",
    )
    fixperms_parser.add_argument(
        "--remote-data-root",
        default=str(DEFAULT_REMOTE_DATA_ROOT),
        help=f"Server-side root used for the rig folder chmod step (default: {DEFAULT_REMOTE_DATA_ROOT})",
    )
    fixperms_parser.add_argument(
        "--permission-user",
        default="",
        help="Optional SSH username for the permission repair host",
    )
    fixperms_parser.add_argument(
        "--permission-host",
        default=DEFAULT_PERMISSION_HOST,
        help=f"SSH host used for permission repair (default: {DEFAULT_PERMISSION_HOST})",
    )

    return parser


def main(argv=None):
    ensure_logging_configured()
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "fixdaq":
        transfer_fixdaq(
            source_root=args.source_root,
            delete_zarr=args.delete_zarr,
        )
    elif args.command == "push":
        transfer_push(
            **_resolve_common_kwargs(args),
            repair_permissions=args.repair_permissions,
        )
    elif args.command == "all":
        transfer_all(
            **_resolve_common_kwargs(args),
            repair_permissions=args.repair_permissions,
            delete_zarr=args.delete_zarr,
        )
    elif args.command == "fixperms":
        transfer_fixperms(
            rig_name=args.rig_name,
            remote_data_root=args.remote_data_root,
            permission_user=args.permission_user,
            permission_host=args.permission_host,
        )
    else:
        parser.error(f"Unknown command: {args.command}")
    return 0
