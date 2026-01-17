# Code for "Transient Distortions of the South Atlantic Anomaly Radiation Environments Driven by Electric Fields"

This repository contains the simulation codes used in the study:
> **Transient Distortions of the South Atlantic Anomaly Radiation Environments Driven by Electric Fields**  
> Authors: Ze-Fan Yin et al.

The codes conduct test-particle simulations in an eccentric dipole magnetic field to reproduce the main results and figures presented in the paper.

---

## 📂 Repository structure
-  `20240409/`, `20240924/`,`20240715/` – Source codes for three SAA distortion events, corresponding to Event I, II, and III in the study
- `Backward_MSS_ED_Single.m` – Main script for test-particle simulation  
- `Traj_particle_bounce_drift_ED.m` – Function to trace particle trajectories  
- `Calc_B_ED.m`, `Calc_E_ED.m` – Functions to calculate the eccentric dipole magnetic and electric fields  
- `PreCalc_A_0409.mat`,`PreCalc_B_0924.mat`, `PreCalc_A_0715.mat` – Pre-processed MSS1 observations (raw data available at [Zenodo](https://doi.org/10.5281/zenodo.18277877))  
- `igrfmex_wrapper_GC.m`, `igrfmex.mexmaca64` – Wrapper functions to link the original geopack source files.
  - *Note: Source codes and build scripts are available in the `GeopackMex/` folder. The MEX compilation approach is adapted from [RibomBalt/geopack_mex](https://github.com/RibomBalt/geopack_mex/).*

## ⚙️ Requirements
- MATLAB R2024b  
- Tested on macOS 15.6 (should also work on other OS with MATLAB installed) 
  - *Note: Windows/Linux users should recompile the source files provided in the `GeopackMex/` folder based on the script `rebuild_mex.sh`.*

## 📈 Reproducing paper figures
- **Trajectories of representative electrons in Figs. 5, 9 & Supplementary Fig. 3** – Run `Backward_MSS_ED_Single.m` directly. Modify `T_loop_range` for three different moments (see code comments).
- **Full distribution in Figs. 4, 9 & Supplementary Fig. 3** – Modify the loop ranges (energy, time, and pitch angle) in `Backward_MSS_ED_Single.m` to obtain energy spectrograms and pitch angle distributions

## 📃 License
This code is released under the MIT License. See the `LICENSE` file for details.  
