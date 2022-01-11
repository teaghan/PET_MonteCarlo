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

patient_dir="../patients/${patient_name}/"

# Create phantom
python ../phantom/generate_phantom.py "${patient_name}"

# Split simulation into multiple jobs
python setup_batch_jobs.py "${patient_name}"

# Run the simulations
python start_simulations.py "${patient_name}"
