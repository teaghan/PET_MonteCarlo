##

# Name of patient corresponding to the config file you have created
# in the PETRecML/patients/configs/ directory
patient_name="$1"

# Check if name was provided
if [ ${#patient_name} -lt 1 ]
then
  echo Please provide a patient name
  exit 1
fi

python run_unfinished_simulations.py "${patient_name}"


#patient_dir="../patients/${patient_name}/"

# Create phantom
#python check_root_files.py "${patient_name}"


#cd "${patient_dir}scripts/"
#num_jobs=$(find . -type f | wc -l)
# Compute Canada account defined in config file
#account=$(awk '/account/{print $3}' "../../configs/${patient_name}.ini")
# Run all jobs simultaneously
#python ../../../compute-canada-goodies/python/queue_cc.py --account "${account}" --todo_dir "../scripts" --done_dir "../done" --output_dir "../stdout" --num_jobs $num_jobs --num_runs 1 --num_gpu 0 --num_cpu 1 --mem 24G --time_limit "0-03:00"
