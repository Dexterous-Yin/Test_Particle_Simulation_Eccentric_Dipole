# Code for "Transient Distortions of the South Atlantic Anomaly Radiation Environments Driven by Electric Fields"

This repository contains the simulation codes used in the study "Transient Distortions of the South Atlantic Anomaly Radiation Environments Driven by Electric Fields". The code conducts the test-particle simulations in the eccentric dipole magnetic field, to reproduce the main results and figures presented in the paper.

---

## Repository structure
- `20240409/` & `20240715/`: Source code for two events
- `Backward_MSS_ED_Single.m`: Main script for test-particle simulation; 
- `Traj_particle_bounce_drift_ED.m`: Main function to trace particles
-  `Calc_B_ED.m` & `Calc_E_ED.m`: Functions to determine the magnetic and electric fields
-  `PreCalc_A_0409.mat` & `PreCalc_A_0715.mat`: MSS1 observations, with the raw data available in https://doi.org/10.5281/zenodo.16925320
-  `igrfmex_wrapper_GC.m` & `igrfmex.mexmaca64`: Functions to link the original geopack source files

## Requirements
- Programming language: MATLAB R2024b
- Environment: Mac OS 15.6

## Reproducing paper figures
- Figures 5 & 9: Directly run `Backward_MSS_ED_Single`
- Figures 4 & 8: Expand the loop range for energy, time, and pitch angles in `Backward_MSS_ED_Single` to obtain the energy spectrogram and pitch angle distributions
