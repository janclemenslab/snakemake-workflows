# merge_splits

Creates one HDF5 file that links to per-split HDF5 inputs via HDF5 external links. This is intended for merging outputs from `split_videos`-style workflows without copying the actual datasets.

## Expected rule shape

```yaml
rule merge_splits:
    input:
        split_tracks
    output:
        "res/{directory}/{directory}_tracks.h5"
    wrapper:
        "file:../snakemake-workflows/wrappers/merge_splits"
```

## Inputs

- `input`: one or more HDF5 files to merge

## Outputs

- `output[0]`: merged HDF5 file containing `roi01`, `roi02`, ... external links

## Snakemake params

- none

## Notes

- the wrapper removes an existing output file before recreating it
- inputs are stored in the output file attributes under `inputs`
- each input becomes an external link named `roiNN`
