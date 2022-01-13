#!/bin/bash

# System specific module loads 
#source $HOME/project/obriaint/medphys/bin/activate
module load StdEnv/2016.4 
module load arch/avx2
module load nixpkgs/16.09 gcc/5.4.0 gate/9.0
source `which geant4.sh`

# Create root files
hadd combined.root 1.root 2.root
