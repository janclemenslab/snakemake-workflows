# Snakemake workflows

## Installation

### User
- you need to have miniforge installed on your machine
- make fab env:
```shell
conda create -n fab python=3.14 uv ffmpeg
conda activate fab
uv pip install git+https://github.com/janclemenslab/snakemake-workflows
```
- setup passwordless ssh for the cluster: [https://hpcwiki.uol.de/hpc-usage/the-basics/login/ssh-keys]

#### Note on local installs (ignore)
- if this checkout lives under a path containing `#` (for example `/Volumes/agauneu/#Data/...`), `uv pip install -e .` currently fails because `uv` misparses the local path; use `pip install -e .` locally or clone the repo to a path without `#`


### Usage
- start local session in the terminal on your PC
- go to the folder for your experiments (chainigmic, backlight etc) - on windows: `cd W:/#Data/chainingmic`
- annotate if necessary: `fab annotate`
- submit jobs: `fab submit --user abcd1234`
- monitor the latest controller run in a compact summary: `fab monitor --user abcd1234`
- open a local FastHTML dashboard for the latest controller run: `fab dashboard --user abcd1234`
- check status of tracker jobs: `fab queue --user abcd1234`

#### Logging, monitoring, and submission

- `fab submit --user abcd1234`
  - runs Snakemake remotely on `rosa.hpc.uni-oldenburg.de` inside the shared cluster env
  - starts the controller in the background with `nohup`
  - writes one controller log per launch under `log/slurm/controller-YYYYMMDDTHHMMSS.log`
  - runs `--conda-create-envs-only` first and then launches the actual workflow
- `fab monitor --user abcd1234`
  - parses the newest controller log under `log/slurm/`
  - shows a compact terminal summary of progress, submitted jobs, failed jobs, and local rules
  - merges in live `squeue` state when the selected run is still active
- `fab dashboard --user abcd1234`
  - starts a local FastHTML server, by default on [http://127.0.0.1:2718/](http://127.0.0.1:2718/)
  - reads local project files such as `log/slurm/controller-*.log`, rule logs, and SLURM logs
  - uses the supplied SSH user/host to query live queue state and to submit or cancel runs
  - opens the browser automatically unless launched in headless mode
- `fab queue --user abcd1234`
  - runs `squeue -u <user>` on the cluster for a raw queue view

#### Dashboard for logging and submission

The dashboard has two views:

- `Logs`
  - choose any controller log under `log/slurm/`
  - inspect summary cards, job tables, state/rule filters, and automatic refresh
  - open a job modal that shows controller context, rule log output, and the matching SLURM log
  - open a workflow DAG modal derived from the selected controller log
  - cancel the selected SLURM run directly from the page when the run is still active
- `Submit`
  - `Submit All` is equivalent to `fab submit` with the workflow default target
  - `Selective Submit` builds explicit Snakemake targets and submits them with `--force`
  - for `playback`, the curated rules are `detect_fly_chambers`, `merge_tracks`, `pb_speed`, and `report_plots`
  - for `chainingmic`, the curated rules are `sleap` and `song`
  - for other workflows, the dashboard falls back to discovered rule names from the Snakefile but does not yet build curated experiment checklists from `dat/`

Useful `fab dashboard` task options:

- `--port 2720` to change the local port
- `--bind-host 0.0.0.0` to serve beyond localhost
- `--remote-host some.cluster.example` to override the default SSH host
- `--limit 200` to change the number of jobs shown per page
- `headless`, `open_browser`, and `verbose` are also supported task parameters

### Transfer from Windows
- stage incoming folders under `~/data.transfer`
- after installing `snakemake-workflows`, use the `transfer` command directly
- add `--permission-user abcd1234` only if the SSH login user should differ from the default on that machine
- run `transfer all --rig-name chainingmic` to convert DAQ zarr folders first, then move and repair the rig folders
- run `transfer push --rig-name chainingmic` to move the folders into the SMB share and then quietly run `chmod g+rwx` and `chmod a+rx` on `<rig>/dat` and `<rig>/res` over SSH on `rosa.hpc.uni-oldenburg.de`
- rerun `transfer fixperms --rig-name chainingmic` to reapply those chmods on the rig's `dat` and `res` folders

### If you want to play with the env on the cluster
login with `ssh abcd1234@rosa.hpc.uni-oldenburg.de`,
then:
```shell
umask a=rwx
module load Miniforge3
conda activate "/fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb"
cd "/fs/s6k/groups/agauneu/#Data"
```



### Admin (ignore me)
Make shared `ncb` environment to be used by everyone on the cluster:
```shell
umask a=rwx
module load Miniforge3
conda create --prefix "/fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb" python=3.14 git uv ffmpeg
conda activate "/fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb"
uv pip install snakemake fabric rich pyyaml snakemake-executor-plugin-slurm h5py matplotlib av
uv pip install -e /fs/s6k/groups/agauneu/#Data/snakemake-workflows
```



## Wrappers
### How to make a new one
See the repo and the official [snakemake-wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/index.html#contribute) for examples.

Each workflow contains:
     - `wrapper.py` - creates the output for that rule. Has access to params, input, outputs etc of the snakefile of the rule (and of all other rules).
    - `environment.yaml` - conda environment for the wrapper. Will be created automatically by snakemake and used for executing the rule. Snakemake will update the environment when the rule file changes.
    - `meta.yaml` (optional) - contains documentation like the expected inputs and outputs, required parameters etc.
    - tf models or other stuff required by the job


### Existing wrappers

Only wrapper directories that contain a `wrapper.py` are runnable Snakemake wrappers. Right now those are:

| Wrapper | Purpose | Typical input | Typical output | Details |
| --- | --- | --- | --- | --- |
| [`sleap`](wrappers/sleap/README.md) | Run SLEAP tracking on a video, trim unreadable tail frames if necessary, and export HDF5 tracks. | `dat/<directory>/<sample>.mp4` | `res/<directory>/<sample>_sleap.h5` | See wrapper README for model and tracking params. |
| [`das`](wrappers/das/README.md) | Run DAS on DAQ recordings and export annotations as CSV. | `dat/<directory>/<directory>_daq.h5` | `res/<directory>/<directory>_annotations.csv` | Supports one model directory or a list of model directories. |
| [`split_videos`](wrappers/split_videos/README.md) | Detect fly chambers, decide which chambers contain a fly, and crop one MP4 per active chamber. | `dat/<directory>/<sample>.mp4` | `dat/<directory>/<sample>_fly_chambers/` | Writes a detections manifest, debug image, and per-chamber videos. |
| [`merge_splits`](wrappers/merge_splits/README.md) | Combine per-split HDF5 outputs using HDF5 external links without copying datasets. | one or more split HDF5 files | merged HDF5 file such as `res/<directory>/<sample>_tracks.h5` | Each input becomes an external link like `roi01`, `roi02`, ... |
| [`pb_speed`](wrappers/pb_speed/README.md) | Convert playback tracks into stimulus-aligned speed traces plus a playlist CSV. | merged SLEAP tracks plus chamber manifest | `*_spd.npz` and `*_playlist.csv` | Used in the playback workflow after chamber detection and tracking. |

`wrappers/pb_tune/` currently contains scratch notebooks and style assets only; it is not a runnable wrapper because it has no `wrapper.py`.


## Analysis profiles
Are basically form definition files for formbuilder.py. Params should start with the workflow folder name.

```yaml
- name: workflow.param1
  label: workflow parameter 1
  type: bool
  default: True
- name: sleap.tracking.match
  label: sleap tracking matching method
  type: list
  default: hungarian
  options: hungarian,greedy
- name: sleap.tracking.track_window
  label: sleap tracking window
  type: int
  default: 15
```

Will create the analysis.yaml file (assuming all parameters were unchanged in the GUI):
```yaml
rule_name:
    param_name: True
sleap:
    tracking.match: hungarian
    tracking.track_window: 5
```

Can be accessed in the wrapper.py via: `params = snakemake.params[0]`. `params[snakemake.rule]` will retrieve the associated rule's parameters (`snakemake.rule` is the rulename in the snakefile).



Will use `workflow/analysis_profiles/default.yaml` if that file does not exist. Note that the default schema needs to have an entry for each rule - otherwise the rule will not be run. Right now you have to create a "dummy" entry that starts with the rulename followed by a `.`. The following will run the `pb_speed` rule. This is not great since this will create a check box in the annotator GUI, giving the illusion that this is configurable, but the param value will be ignored and `pb_speed` will run irrespective of the state of this field.
```yaml
- name: pb_speed.
  label: postprocess tracks (dummy option - will be run regardless of selection)
  type: bool
  default: True
```
