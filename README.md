# Snakemake workflows
For use as wrappers in snakefiles


### How to make a new one
See the repo and the official [snakemake-wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/index.html#contribute) for examples.

Each workflow contains:
     - `wrapper.py` which creates the output for that rule from the inputs
    - `environment.yaml` - conda environment for the wrapper. Will be created automatically by snakemake and used for executing the rule. Snakemake will update the environment when the rule file changes.
    - `meta.yaml` (optional)
    - tf models or other stuff required by the job
