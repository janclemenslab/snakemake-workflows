# Snakemake workflows
For use as wrappers in snakefiles


### How to make a new one
See the repo and the official [snakemake-wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/index.html#contribute) for examples.

Each workflow contains:
     - `wrapper.py` which creates the output for that rule from the inputs
    - `environment.yaml` - conda environment for the wrapper. Will be created automatically by snakemake and used for executing the rule. Snakemake will update the environment when the rule file changes.
    - `meta.yaml` (optional)
    - tf models or other stuff required by the job


### Troubleshooting bugs in snakemake and conda
- For some envs (sleap), snakemake can't identify the python because (in the case of sleap) `python --version` does not return in the expected format. Just override this in `snakemake/script.py:
```python
 def _get_python_version(self):
        out = self._execute_cmd("python --version", read=True).strip()
        if not len(out):
             out = "Python 3.5"
        return tuple(map(int, PY_VER_RE.match(out).group("ver_min").split(".")))
```
- I had trouble testing this on windows, with the conda envs being created somewhere else - this is because of the way the env creation command is assembled by snakemake (in snakemake/deployment/conda.py:line304) - the path is quoted, which conda interprets as a relative path. Remove the quotes in  conda.py:line304 - should like like `"--file {}".format...` and `"--prefix {}".format...`.