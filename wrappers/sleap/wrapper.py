import os
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


def resolve_legacy_model_dir(model_dir: str) -> str:
    legacy_prefix = "../snakemake-workflows/sleap/models"
    new_prefix = "../snakemake-workflows/wrappers/sleap/models"
    if model_dir.startswith(legacy_prefix):
        candidate = model_dir.replace(legacy_prefix, new_prefix, 1)
        if os.path.exists(candidate):
            return candidate
    return model_dir


def main():

    from sleap_io.io.main import load_file

    print(f"Testing {snakemake.input.video}")  # noqa: F821
    mv = load_file(snakemake.input.video)  # noqa: F821
    nb_frames = mv.shape[0]

    print(f"  Video has {nb_frames} frames.")

    failed_frames = []
    max_back = 2_000
    print(f"   Testing frames {nb_frames - max_back}-{nb_frames}.")
    for back in range(max_back):
        try:
            mv[nb_frames - back]
        except Exception:
            print(f"        Failed to read frame {nb_frames - back}.")
            failed_frames.append(nb_frames - back)

    nb_frames_corrected = min(failed_frames) - 1
    correction = nb_frames - nb_frames_corrected

    # Monkey-patch SLEAP video IO to ignore unreadable tail frames.
    # @property
    # def num_frames(self, value):
    #     self.num_frames = value

    # @property
    # def last_frame_idx(self):
    #     return int(self.frames - correction - 1)
    # mv.backend.__class__.num_frames = nb_frames - correction - 1
    # sleap.io.video.MediaVideo.last_frame_idx = last_frame_idx

    print(
        f"   Could read only up to {nb_frames_corrected} frames. Will ignore the rest."
    )

    params = snakemake.params[0][snakemake.rule].copy()  # noqa: F821

    print(params)
    # modelname refers to the dir that contains the folder of the individual models
    # for a sleap pipeline (e.g. centroid and centered instance models)
    model_dir = resolve_legacy_model_dir(params.pop("modelname"))
    models = []
    for d in os.listdir(model_dir):
        model_path = os.path.join(model_dir, d)
        if os.path.isdir(model_path):
            models.append(model_path)

    extra = [f"--model_paths {model} " for model in models]

    # normalize params from old analysis files - TODO: update workflow spec and old analysis files
    print("Normalizing params")
    for key, value in params.items():
        if key == "tracking.pre_cull_to_target":
            value = 1 if value else 0
        if key == "no-empty-frames":
            if not value:
                continue
            value = ""
        elif key == "tracking.similarity":
            continue

        if key == "tracking.match":
            key = "track_matching_method"
        elif key == "tracking.robust":
            key = "robust_best_instance"
        # elif key == "tracking.similarity":
        #     key = "scoring_method"
        elif key == "tracking.track_window":
            key = "tracking_window_size"
        elif key == "tracking.tracker":
            if value == "simple":
                continue
            else:
                key = "use_flow"
                value = ""

        key = key.replace(".", "_")
        key = key.replace("-", "_")

        extra.append(f"--{key} {value} ")

    extra = "".join(extra)

    print(extra)

    for out in snakemake.output:
        print(f"Creating {out}")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        cmd = f"export PYTHONIOENCODING=utf-8; sleap track --data_path {snakemake.input.video} --frames 0-{nb_frames_corrected} --output_path {out}.slp --tracking {extra}"
        print(cmd)
        os.system(cmd)
        cmd = f"export PYTHONIOENCODING=utf-8; sleap export {out}.slp -o {out}"
        print(cmd)
        os.system(cmd)


logfile = snakemake.log[0] if snakemake.log else None  # noqa: F821
if logfile:
    log_path = Path(logfile)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as handle, redirect_stdout(handle), redirect_stderr(handle):
        main()
else:
    main()
