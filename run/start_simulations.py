import numpy as np
import os
import sys
import configparser

# Command line arguments
patient_name = sys.argv[1]

# Patient configuration
cur_dir = os.path.dirname(os.path.realpath(__file__))
config_dir = os.path.join(cur_dir, '../patients/configs/')
config = configparser.ConfigParser()
config.read(os.path.join(config_dir, patient_name+'.ini'))

# Number of separate jobs to launch
if len(sys.argv)>2:
    num_jobs = int(sys.argv[2])
else:
    num_jobs = int(config['Simulation']['num_jobs'])
# Compute Canada account to run jobs under
account = config['Simulation']['account']
# Time per job
job_time = config['Simulation']['job_time']
# Memory per job in GB
mem = config['Simulation']['mem']

patient_dir = os.path.join(cur_dir, '../patients/', patient_name)
data_dir = os.path.join(patient_dir, 'phantom')
results_dir = os.path.join(patient_dir, 'results')
scripts_dir = os.path.join(patient_dir, 'scripts')

# Compute-canada goodies command
cmd = 'python ../../../compute-canada-goodies/python/queue_cc.py '
cmd += '--account "%s" --todo_dir "../scripts" ' % (account)
cmd += '--done_dir "../done" --output_dir "../stdout" '
cmd += '--num_jobs %i --num_runs 1 --num_gpu 0 ' % (num_jobs)
cmd += '--num_cpu 1 --mem %sG --time_limit "%s"' % (mem, job_time)

# Make root combination script executable
os.system('chmod u+x %s' % os.path.join(results_dir, 'combine_root.sh'))

# Execute jobs
os.chdir(scripts_dir)
os.system(cmd)
