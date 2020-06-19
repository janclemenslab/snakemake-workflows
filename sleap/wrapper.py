from snakemake.shell import shell
import os

params = snakemake.params[0][snakemake.rule].copy()

# modelname refers to the dir that contains the folder of the individual models
# for a sleap pipeline (e.g. centroid and centered instance models
model_dir = params.pop('modelname')
models = []
for d in os.listdir(model_dir):
    model_path = os.path.join(model_dir, d)
    if os.path.isdir(model_path):
        models.append(model_path)
extra = [f'--model {m} ' for m in models]
extra.extend([f'--{k} {v} ' for k, v in params.items()])

for out in snakemake.output:
    shell("umask g+rwx; \
           export CUDA_VISIBLE_DEVICES=0; \
           sleap-track {snakemake.input.video} --output {out} {extra};")
    # load results and do what sleap-convert does while preserving confidence scores
    # shell("umask g+rwx; \
    #        export CUDA_VISIBLE_DEVICES=0; \
    #        sleap-track {snakemake.input.video} --output {out} {extra}; \
    #        sleap-convert --output {out} --video {snakemake.input.video} --format analysis {out}")
