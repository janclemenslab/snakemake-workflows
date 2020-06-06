from snakemake.shell import shell

params = snakemake.params[0][snakemake.rule].copy()

extra = [f'--model {m} ' for m in params.pop('modelname')]
extra.extend([f'--{k} {v} ' for k,v in params.items()])

for out in snakemake.output:
    shell("umask g+rwx; \
           sleap-track {snakemake.input.video} \
           -output {out} {extra}")