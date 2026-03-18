# split_videos

Detects fly chambers in a multi-chamber video from the first frame, identifies which chambers contain a fly, writes a detections manifest and debug image, and crops one MP4 per active chamber.

## Expected rule shape

```yaml
checkpoint split_videos:
    input:
        video="dat/{directory}/{directory}.mp4"
    output:
        directory("dat/{directory}/{directory}_fly_chambers")
    wrapper:
        "file:../snakemake-workflows/split_videos"
```

## Inputs

- `input.video`: source video (`.mp4`)

## Outputs

- `output[0]`: output directory, usually `dat/{directory}/{directory}_fly_chambers`
- `{stem}_detections.json`
- `{stem}_detections.png`
- `{stem}_chamberNN.mp4` for each chamber with a detected fly

## Snakemake params

- `min_chambers`: minimum number of chambers that must be detected before the wrapper fails. Default `6`.
- `crop_padding`: padding in pixels added around each cropped chamber video. Default `10`.
- `detection_padding`: padding in pixels added around each chamber when searching for the fly. Default `18`.
- `force`: if truthy, overwrite existing cropped videos, manifest, and debug image. Default `False`.
- `skip_cropping`: if truthy, only write the detections manifest and debug image. Default `False`.

## Notes

- the output directory must be next to the input video
- the output directory name must start with the input stem, because the wrapper infers the suffix from it
- the manifest contains one record per detected chamber and a `has_fly` flag that downstream rules can use to expand jobs dynamically
