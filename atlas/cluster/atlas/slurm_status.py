#!/usr/bin/env python3
"""
Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --cluster slurm_scheduler.py --cluster-status slurm_status.py
"""

import os
import sys
import warnings
import subprocess


jobid = sys.argv[1]


out= subprocess.run(['scontrol','show','jobid',jobid],stdout=subprocess.PIPE).stdout.decode('utf-8')

def parse_key_value(stream):
    params={}
    for key_value_pair in stream.split():
        name, var = key_value_pair.partition("=")[::2]
        params[name.strip()] = var
    return params


state=parse_key_value(out)['JobState']


map_state={"PENDING":'running',
           "RUNNING":'running', 
           "SUSPENDED":'running', 
           "CANCELLED":'failed', 
           "COMPLETING":'running', 
           "COMPLETED":'success', 
           "CONFIGURING":'running', 
           "FAILED":'failed',
           "TIMEOUT":'failed',
           "PREEMPTED":'failed',
           "NODE_FAIL":'failed',
           "REVOKED":'failed',
           "SPECIAL_EXIT":'failed',
           "":'success',
	   "OUT_OF_MEMORY":'failed'}

print(map_state.get(state,'failed'))


