from snakemake.shell import shell

params = snakemake.params[0][snakemake.rule]
extra = ''.join([f'--{k} {v} ' for k,v in params.items()])

for out in snakemake.output:
    shell("umask g+rwx; \
           python3 -m ethotracker.tracker.pursuit {snakemake.input.video} {out} {extra}")