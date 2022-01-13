#!/bin/bash

# System specific module loads 
module load StdEnv/2016.4
module load arch/avx2
module load nixpkgs/16.09 gcc/5.4.0 gate/9.0
source `which geant4.sh`

# Run Monte-Carlo simulation
cd /path/to/gate/
Gate -a "[scale,0.0001] [root_fn,path/to/output/root] [mat_fn,AttenuationRange.dat] [atn_fn,xcat_atn.h33] [act_fn,xcat_act.h33] [TimeSlice,1] [TimeStart,0] [TimeStop,10]" main_phantom.mac > output.txt &