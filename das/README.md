# das

Annotates song or vibration recordings with DAS and writes the results as a CSV table.

## Expected rule shape

```yaml
rule das:
    input:
        recording="dat/{directory}/{directory}_daq.h5"
    output:
        "res/{directory}/{directory}_annotations.csv"
    params:
        read_annotations
    wrapper:
        "file:../snakemake-workflows/das"
```

## Inputs

- `input[0]`: recording file used to assemble the experiment dataset

## Outputs

- one CSV file per output entry in `snakemake.output`

## Snakemake params

- `modelname`: path or list of paths to DAS model directories

## Notes

- the wrapper reads rule-specific params from `snakemake.params[0][snakemake.rule]`
- the wrapper loads the recording through `xarray_behave`, runs DAS prediction, converts predictions to timestamps, and saves a CSV
- if `modelname` is a list, predictions from all listed models are merged before export
