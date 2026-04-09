# pb_speed

Converts merged playback SLEAP tracks into stimulus-aligned speed traces saved as `*_spd.npz` and an event playlist saved as `*_playlist.csv`.

## Expected rule shape

```yaml
rule pb_speed:
    input:
        tracks="res/{session}/{video}_sleap.h5"
    output:
        speed="res/{session}/{video}_spd.npz"
        playlist="res/{session}/{video}_playlist.csv"
    log:
        "log/{session}/{video}_pb_speed.log"
    wrapper:
        "file:../snakemake-workflows/wrappers/pb_speed"
```

## Notes

- the wrapper infers `root`, `session`, and `video` directly from the merged SLEAP filename
- DAQ metadata is loaded from `dat/{session}/{session}_daq.log`
- traces are aligned to song-driven stimulus onsets and flattened to `(trials * flies, time)`
- trace rows are ordered in fly-major blocks: all trials for fly 1, then all trials for fly 2, etc.
- the playlist CSV includes the DAQ metadata aligned to detected events, plus `onset_times` and `offset_times`, duplicated to the same `(trials * flies)` row order
