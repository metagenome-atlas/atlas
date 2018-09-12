#!/usr/bin/env python3


import sys
from subprocess import Popen, PIPE

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# let snakemake read job_properties
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)





#default paramters defined in cluster_spec (accessed via snakemake read_job_properties)
cluster_param= job_properties["cluster"]

cluster_param['name'] = job_properties['rule']

# overwrite default parameters if defined in rule (or config file)
if 'threads' in job_properties:
    cluster_param["threads"] = job_properties["threads"]
for res in ['time','mem']:
    if res in job_properties["resources"]:
        cluster_param[res] = job_properties["resources"][res]

#print parameters
#eprint("Submit job with parameters:\n"+"\n".join(["\t{} : {}".format(key,cluster_param[key]) for key in cluster_param]))

# check which system you are on and load command command_options

command_options= cluster_param['command_options'][cluster_param['system']]
command= command_options['command']

key_mapping= command_options['key_mapping']

# construct command:
for  key in key_mapping:

    if key in cluster_param:
        command+=" "
        command+=key_mapping[key].format(cluster_param[key])

command+=' {}'.format(jobscript)

eprint("submit command: "+command)

p = Popen(command.split(' '), stdout=PIPE, stderr=PIPE)
output, error = p.communicate()
if p.returncode != 0:
   raise Exception("Job can't be submitted\n"+output.decode("utf-8")+error.decode("utf-8"))
else:
	print(output.decode("utf-8"))
