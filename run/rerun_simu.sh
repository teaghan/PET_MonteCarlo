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