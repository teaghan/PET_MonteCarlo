# PET Monte Carlo
Simulation and analysis code to accompany the paper (still to come...)

This repository includes scripts to perform Monte Carlo simulations of the Siemens Biograph 40 mCT PET/CT system using the GEANT4 application for tomographic emission (GATE).

## Dependencies

- [GATE](https://opengate.readthedocs.io/en/latest/installation.html)
  - this includes install ROOT and GEANT4
- SciPy: `pip install scipy`
- scikit-image: `pip install scikit-image`

## Running the GATE simulations

1. Before running a simulation, you may have to adjust the module loads for a few of the scripts ([run/mc_job.sh](./run/mc_job.sh) and [run/combine_root.sh](./run/combine_root.sh)) so that they properly load the correct modules on your system. If you have GATE and so forth already loaded then you can just remove these lines from the top of the scripts.
2. For each simulation, you will need to create a **configuration file** in the [config directory](./patients/configs/). For example, take a look at this [config file](./patients/configs/line_10.ini) for a line source simulation. You will notice that there is a parameter set as `num_jobs=100`. This means that the simulation is going to be split into 100 separate jobs.
3. Using this configuration, you can generate the voxelized phantom using the commmand `python phantom/generate_phantom.py line_10` which will create a directory tree for this patient in the [patient directory](./patients/).
4. Next, you can create scripts for each job using the command `python setup_batch_jobs.py line_10` which will save the scripts within the newly created `patients/line_10/scripts/` directory.
5. If you are running the scripts on one of the Compute Canada servers, you can run the command `python start_simulations.py line_10` to start these simulations, otherwise, you will have to run these separate jobs using your own method.

All of these steps can be accomplished by running the ([run/run_simu.sh](./run/run_simu.sh) script using the command `./run_simu.sh line_10` inside the [run directory](./run/).

## Analysis code

The code required to bin the root output into sinograms and reconstruct the images could not be made publicly available, so there is a slight jump from the GATE simulation outputs to the sinograms and images. Therefore, we have provided the necessary data files for download.

### Data download

The data can be downloaded [here](https://zenodo.org/record/5851646) or possibly by just clicking this [link](https://zenodo.org/record/5851646/files/validation_data.zip?download=1).

Once downloaded, unzip the file and place each of the subdirectories within the [patient directory](./patients/). For instance, after doing this you should have a directory tree `PET_MonteCarlo/patients/nema_simu`.

### Analysis notebooks

Checkout the following notebooks to see how we compared the measured data to our simulations:

- [Uniform cylindrical phantom](./analysis/uniform_cylinder.ipynb)
- [Line sources](./analysis/line_source.ipynb)
- [NEMA IEC body image quality phantom](./analysis/nema_image_quality.ipynb)

