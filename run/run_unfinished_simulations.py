import os
import sys
import uproot
import configparser
import glob

# Command line arguments
patient_name = sys.argv[1]

# Directories
cur_dir = os.path.dirname(os.path.realpath(__file__))
patient_dir = os.path.join(cur_dir, '../patients/', patient_name)
results_dir = os.path.join(patient_dir, 'results')
config_dir = os.path.join(cur_dir, '../patients/configs/')
results_dir = os.path.join(patient_dir, 'results')
scripts_dir = os.path.join(patient_dir, 'scripts')

# Patient configuration
config = configparser.ConfigParser()
config.read(os.path.join(config_dir, patient_name+'.ini'))

# Number of separate jobs to launch
num_jobs = int(config['Simulation']['num_jobs'])
# Compute Canada account to run jobs under
account = config['Simulation']['account']
# Time per job
job_time = config['Simulation']['job_time']
# Memory per job in GB
mem = config['Simulation']['mem']

#'''
print('Checking %i jobs' % (num_jobs))
for n in range(1, num_jobs+1):
    # Determine if the current job needs to be rerun
    rerun = False
    delete_root = False
    root_fn = os.path.join(results_dir, '%i.root' % n)
    if os.path.isfile(root_fn):
        # Read current sinogram
        try:
            root = uproot.open(root_fn)
            if ((len(root.keys())>8) or (len(root.keys())==0)):
                print('Rerunning job %i' % (n))
                rerun = True
                delete_root = True
        except:
            rerun = True
            delete_root = True

    else:
        print('Root file for job %i does not exist. Rerunning...' % (n))
        rerun = True
        
    if rerun is True:
        script_fn = glob.glob(os.path.join(patient_dir, 'done/', '*_%i.sh' % (n)))
        if len(script_fn)>0:
            os.system('mv %s %s' % (script_fn[0],
                                    os.path.join(patient_dir, 'scripts/')))
    if delete_root is True:
        os.system('rm %s' % (root_fn))
    if (n+1)%500==0:
        print('%i complete' % (n+1))
#'''
        
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
