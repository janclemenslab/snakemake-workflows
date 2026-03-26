# sleap

Runs SLEAP tracking on a video, trims unreadable tail frames if needed, and exports the final result as HDF5.

## Expected rule shape

```yaml
rule sleap:
    input:
        video="dat/{directory}/{directory}.mp4"
    output:
        "res/{directory}/{directory}_sleap.h5"
    params:
        read_annotations
    log:
        "log/{directory}/{directory}_sleap.log"
    wrapper:
        "file:../snakemake-workflows/wrappers/sleap"
```

## Inputs

- `input.video`: source video (`.mp4`)

## Outputs

- one exported HDF5 file per output entry in `snakemake.output`
- an intermediate `.slp` file is also created next to each output before export

## Snakemake params

The wrapper reads params from `snakemake.params[0][snakemake.rule]`.

Common params from `meta.yaml`:
- `modelname`: path to a directory containing one or more SLEAP model subdirectories
- `tracking.tracker`
- `tracking.clean_instance_count`
- `tracking.similarity`
- `tracking.match`
- `tracking.track_window`

Additional SLEAP CLI options can be passed through the same param block. The wrapper normalizes several legacy keys before building the command line.

## Notes

- the wrapper probes the tail of the video and trims unreadable final frames from the tracking range
- every model subdirectory inside `modelname` is passed as a `--model_paths` argument
- if a Snakemake log file is provided, stdout and stderr are redirected there
