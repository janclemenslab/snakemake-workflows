# AGENTS.md

This repository contains shared Snakemake workflow code, wrapper rules, and Fabric commands used from experiment/project folders under `/Volumes/agauneu/#Data`.

## Keep this file current

- If you discover repository-specific information that is likely to matter for future tasks, update this `AGENTS.md` as part of the same job.
- Record durable facts such as:
  - changes in how code should be run, debugged, or validated
  - important code structure, entry points, path assumptions, or environment constraints
  - workflow conventions that are easy to miss but affect correct edits
- Prefer updating or replacing outdated guidance instead of adding contradictory notes.
- Do not add transient task notes, one-off debugging observations, or temporary workarounds unless they have become a stable part of the workflow.

## Preferred environment

- Use the local conda environment `fab` for Codex work:
  - `conda run -n fab ...`
  - or `conda activate fab`
- Use `fab` for Python scripts, `pyav`, `ffmpeg`, quick Snakemake/debug runs, and wrapper validation.
- Do not rely on `.venv` or wrapper-specific conda envs for ad hoc debugging unless the task explicitly requires it.
- Remote cluster jobs use the shared env at `/fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb`; that is the runtime env for `fab submit`, not the default local dev env.

## Repository structure

- `README.md`
  - user/admin setup notes and high-level usage
- `profiles/rosa/config.yaml`
  - Snakemake workflow profile used for cluster submission on `rosa`
- `src/fab_commands/`
  - Fabric tasks exposed by the package
  - includes `submit`, `monitor`, `dashboard`, `queue`, `unlock`, `fixvideos`, `fixdaq`
  - `monitor.py` parses controller logs and renders terminal summaries
  - `monitor_dashboard.py` is the local FastHTML app behind `fab dashboard`
  - `submission.py` builds workflow-aware dashboard submission targets
- `wrappers/`
  - shared rule wrappers
  - current wrapper folders: `das`, `merge_splits`, `pb_speed`, `sleap`, `split_videos`
  - each wrapper normally contains:
    - `wrapper.py`
    - `environment.yaml`
    - `meta.yaml`
    - optional `README.md`, models, helper files, notebooks
- `.envs/`
  - cached/shared environments; do not edit unless explicitly asked

## Workflow model

- This repo is usually not the working directory for actual Snakemake runs.
- Real workflow runs happen from a project folder such as `/Volumes/agauneu/#Data/<project>`.
- From the project folder, users run Fabric commands like:
  - `fab submit --user <cluster_user>`
  - `fab monitor --user <cluster_user>`
  - `fab dashboard --user <cluster_user>`
  - `fab queue --user <cluster_user>`
- `fab submit` launches the remote Snakemake controller in the background and writes controller logs to `log/slurm/controller-YYYYMMDDTHHMMSS.log`.
- `fab monitor` parses the newest controller log and merges live `squeue` state when the run is still active.
- `fab dashboard` starts a local FastHTML server with:
  - a logs page for controller/rule/SLURM log inspection, queue-aware job tables, and DAG viewing
  - a submit page for workflow-aware remote submission
- Dashboard selective submit is curated for:
  - `playback`: `detect_fly_chambers`, `merge_tracks`, `pb_speed`, `report_plots`
  - `chainingmic`: `sleap`, `song`
- Other workflows fall back to generic rule discovery from the Snakefile and do not yet expose curated experiment checklists in the dashboard.
- Remote submit uses:
  - project root on cluster: `/fs/s6k/groups/agauneu/#Data/<project>/`
  - workflow profile: `../snakemake-workflows/profiles/rosa`
  - shared cluster env: `/fs/s6k/groups/agauneu/#Data/snakemake-workflows/.envs/ncb`
  - the project-local `Snakefile`; `fab submit` does not pass `--snakefile`

## Important path assumptions

- Snakemake execution happens from the project folder, not this repo root.
- Project folders provide the `Snakefile` and project-local paths like:
  - `dat/<recording>/...`
  - `res/<recording>/...`
  - `log/slurm/...`
  - `workflow/analysis_profiles/default.yaml`
  - `workflow/analysis_profiles/default_sleap_boxes.yaml`
- Those `workflow/analysis_profiles/...` files are outside this repo and belong to the project/workflow data directory.
- `playback/workflow/snakefile` is self-contained for YAML/default-parameter loading and no longer depends on `../../snakemake-workflows/scripts/main.smk`.

## Wrapper conventions

- `wrapper.py` files typically end with an unconditional `main()` call.
- Do not import wrapper modules directly for inspection unless you want them to execute immediately.
- For local analysis/debugging, load helper functions by reading the file and stripping the trailing `main()` call, or run the wrapper through Snakemake-compatible code paths.
- When changing wrapper behavior, keep `environment.yaml` and `meta.yaml` aligned with the wrapper inputs/outputs/params.
- Prefer validating wrappers against real files under `/Volumes/agauneu/#Data/...`.
- If a project Snakefile cannot be evaluated locally because a shared include is missing, validate the wrapper with a minimal one-rule Snakefile that points directly at the wrapper and uses concrete absolute input/output paths.
- For wrapper smoke tests, write outputs under a project `debug/` directory so local validation does not overwrite canonical `res/` artifacts.

## Local development notes

- This checkout lives under a path containing `#`, so `uv pip install -e .` can misparse the path.
- Prefer `pip install -e .` locally if an editable install is needed from this checkout.
- For quick repo inspection:
  - use `rg` for code/config search
  - use `conda run -n fab python ...` for Python-based validation
- `conda run -n fab snakemake ...` works for local dry-runs and wrapper smoke tests.
- A reliable local wrapper validation pattern is:
  - create a minimal temporary Snakefile for the wrapper
  - run `conda run -n fab snakemake -s /absolute/path/to/test.smk -n --cores 1`
  - then run the same command without `-n` against isolated debug outputs

## Coding guidance for this repo

- Treat `wrappers/*/scratch*.ipynb` as exploratory unless explicitly asked to edit them.
- Avoid editing generated env folders, cache folders, or data outputs unless the task requires it.
- Be careful with long `ffmpeg` writes on the network mount under `/Volumes/agauneu`; writing to a local temp file and moving it into place is safer.
- When changing `src/fab_commands/__init__.py`, keep the command names and remote path assumptions consistent with `README.md` and `profiles/rosa/config.yaml`.
