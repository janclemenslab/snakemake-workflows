# Snakemake workflows

## Installation

### User
- you need to have miniforge installed on your machine
- make fab env:
```shell
conda create -n fab python=3.14 uv ffmpeg
conda activate fab
uv pip install snakemake fabric rich pyyaml uv snakemake-executor-plugin-slurm git+https://github.com/janclemenslab/etho_annotator

```
- setup passwordless ssh for the cluster: [https://hpcwiki.uol.de/hpc-usage/the-basics/login/ssh-keys]


### Usage
- start local session in the terminal on your PC
- go to the folder for your experiments (chainigmic, backlight etc) - on windows: `cd W:/#Data/chainingmic`
- annotate: `fab annotate`
- convert vids, submit jobs: `fab fixvideos submit -u USERNAME`
- check status of tracker jobs: `fab queue -u USERNAME`

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
uv pip install snakemake fabric rich pyyaml snakemake-executor-plugin-slurm
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

#### SLEAP
[sleap README](sleap/README.md)

#### DAS
[das README](das/README.md)

#### split_videos
[split_videos README](split_videos/README.md)

#### merge_splits
[merge_splits README](merge_splits/README.md)



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



Will use `workflow/analysis_profiles/default.yaml` if that file does not exist. Note that the default schema needs to have an entry for each rule - otherwise the rule will not be run. Right now you have to create a "dummy" entry that starts with the rulename followed by a `.`. The following will run the `spd` rule is run. This is not great since this will create a check box in the annotator GUI, giving the illusion that this is - but the param value will be ignored and the `spd` irrespective of the state of this field.
```yaml
- name: spd.
  label: postprocess tracks (dummy option - will be run regardless of selection)
  type: bool
  default: True
```
