from snakemake.shell import shell

params = snakemake.params[0][snakemake.rule].copy()

extra = [f'--model {m} ' for m in params.pop('modelname')]
extra.extend([f'--{k} {v} ' for k,v in params.items()])

for out in snakemake.output:
    shell("umask g+rwx; \
           module load cuda10.0/toolkit/10.0.130; \
           module load cuda10.0/blas/10.0.130; \
           module load cudnn/10.0v7.6.3; \
           export CUDA_VISIBLE_DEVICES=0; \
           sleap-track {snakemake.input.video} -output {out} {extra}")