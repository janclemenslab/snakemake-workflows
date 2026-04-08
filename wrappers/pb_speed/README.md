# pb_speed

Converts merged playback SLEAP tracks into stimulus-aligned speed traces saved as `*_spd.npz`.

## Expected rule shape

```yaml
rule pb_speed:
    input:
        tracks="res/{session}/{video}_sleap.h5"
    output:
        speed="res/{session}/{video}_spd.npz"
    log:
        "log/{session}/{video}_pb_speed.log"
    wrapper:
        "file:../snakemake-workflows/wrappers/pb_speed"
```

## Notes

- the wrapper infers `root`, `session`, and `video` directly from the merged SLEAP filename
- DAQ metadata is loaded from `dat/{session}/{session}_daq.log`
- traces are aligned to song-driven stimulus onsets detected from the assembled dataset
