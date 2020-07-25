"""Global fab commands."""
from fabric.api import *
import os
import glob


# define globals
NAME: str
REMOTE: str
SERVER: str
FOLDER: str
SCRATCH: str

LOCAL: str = './'
EXCLUDES_PUSH = ["dat/**/*.avi", "dat/**/*.h264"]
EXCLUDES_PULL = []


def update_globals():
    """Call in fabfile to register all funs."""
    global SERVER, REMOTE, FOLDER, SCRATCH
    FOLDER = f"/scratch/clemens10/{NAME}/"
    SERVER = f"{env.user}@{env.hosts[0]}:"
    REMOTE = SERVER + FOLDER
    SCRATCH = "scratch" if "102" in env.hosts[0] else "scratch2"


def _rsync(source, target, delete=False, excludes=[], fix_perms=False):
    RSYNC = 'rsync -rvltHhz --update --progress'
    if delete:  # delete files on target
        RSYNC += ' --delete'

    if fix_perms:  # set perms to read and write (and execute?) for all
        RSYNC += ' --no-p --no-g --chmod=a=rwX'

    for exclude in excludes:  # exclude files/directories
        RSYNC += ' --exclude="{0}"'.format(exclude)

    RSYNC += " {0} {1}".format(source, target)
    print(source, target)
    local(RSYNC)
    if fix_perms:
        with settings(hide('warnings'), warn_only=True):
            run(f"umask a+r,a+w,a+x; chmod ugo=rwX -Rf {FOLDER}{source}/")


def push(delete=False, excludes=EXCLUDES_PUSH, fix_perms=True):
    with settings(warn_only=True):
        _rsync(LOCAL + 'workflow', REMOTE, delete, excludes, fix_perms)
        _rsync(LOCAL + 'dat', REMOTE, delete, excludes, fix_perms)
        _rsync('../snakemake-workflows', SERVER + '/scratch/clemens10/', delete=delete, fix_perms=fix_perms)


def pull(delete=False, excludes=EXCLUDES_PULL):
    _rsync(REMOTE + 'log', LOCAL, delete, excludes)
    _rsync(REMOTE + 'res', LOCAL, delete, excludes)


def _fixvideo(pathname, overwrite=False):
    trunk = os.path.splitext(pathname)[0]
    with settings(warn_only=True):
        if not overwrite and not os.path.exists(f'{trunk}.mp4'):
            local(f'ffmpeg -i "{trunk}.avi" -c copy {trunk}.mp4')  # && rm -vf {trunk}.avi')


def fixvideos(pathname_pattern="dat/**/*.avi", overwrite=False):
    dirnames = glob.glob(pathname_pattern)
    print(dirnames)
    for dirname in dirnames:
        try:
            _fixvideo(dirname, overwrite)
        except Exception as e:
            print(e)


def status():
    run("sacct --format='JobID,Elapsed,MaxVMSize,MaxRSS,CPUTime,State'")


def queue():
    run(f"squeue -u {env.user}")


def submit():
    run(f"umask a=rwx; \
        cd {FOLDER};\
        snakemake --cores 1 --rerun-incomplete --use-conda --conda-create-envs-only --conda-prefix ~/snake_envs;\
        sbatch -t 48:00:00 -o ./log/slurm-%j.out -C {SCRATCH} ../snakemake-workflows/scripts/bsubmit_block.sh {NAME}")