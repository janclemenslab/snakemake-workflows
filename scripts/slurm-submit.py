#!/usr/bin/env python3
"""Submits snakefile jobs via slurm.

How to call:
  python slurm-submit.py [DEPENDENCIES] DATADIR LOGDIR JOBSCRIPT

"""
import os
import sys
import re
import datetime
import subprocess
from snakemake.utils import read_job_properties
import logging


# parse cmd line args
DEPENDENCIES = sys.argv[1:-3]
DATADIR = sys.argv[-3]
LOGDIR = sys.argv[-2]
JOBSCRIPT = sys.argv[-1]
# get all props from snakefile
props = read_job_properties(JOBSCRIPT)

# setup up logging to file
logging.basicConfig(level=logging.DEBUG, filename=f'{LOGDIR}/slurm-submit.log', filemode='a')
logging.debug("COMMAND LINE ARGS:")
logging.debug(sys.argv)

# get info about job from sys
mo = re.match(r'(\S+)/snakejob\.\S+\.(\d+)\.sh', JOBSCRIPT)
assert mo
sm_tmpdir, sm_jobid = mo.groups()

# set up rule and job name
rule = props["rule"]
jobname = f"{rule}-{sm_jobid}"
cmdline = f'sbatch -J {jobname} '

defaults = {'partition': 'medium', 'nodes': "1", 'ntasks': "1", 'time': "48:00:00", 'mem': None,
            'constraint': 'scratch', 'qos': None, 'gpus': None}
defaults.update((k, props["params"][k]) for k in defaults.keys() & props["params"].keys())  # update existing keys with values from snakefile
for key, val in defaults.items():
    if val is not None:
        cmdline += f'--{key} {val} '

# log into LOGDIR with default name or from snakefile params
if props["params"].get("output"):
    logfilename = f'{LOGDIR}/{props["params"].get("output")}'
else:
    logfilename = f"{LOGDIR}/{rule}_{jobname}_output.txt"
cmdline += f"-o {logfilename} "


# attach job dependencies
if DEPENDENCIES:
    cmdline += " -d afterok:{} ".format(":".join(DEPENDENCIES))

# the actual script to run
cmdline += JOBSCRIPT

# the success file
cmdline += f" {sm_tmpdir}/{sm_jobid}.jobfinished "
logging.debug("COMMAND:")
logging.debug(cmdline)

# submit the job and capture output with jobid
logging.debug('CALLING COMMAND:')
try:
    res = subprocess.run(cmdline, check=True, shell=True, stdout=subprocess.PIPE)
    logging.debug(res)
except subprocess.CalledProcessError as e:
    logging.error(e)
    raise e

# parse jobid from sbatch output - print jobid
res = res.stdout.decode()
try:
    m = re.search("Submitted batch job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
    logging.debug(f'   JOBID: {jobid}')
except Exception as e:
    logging.error(e)
