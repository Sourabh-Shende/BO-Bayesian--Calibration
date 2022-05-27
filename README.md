# BO-Bayesian Calibration  
Bayesian optimization (BO) code for the Bayesian calibration of hyperelastic material models of porcine meniscus.

## Installation
This code requires following softwares:
* MATLAB (Numerical computing platform)
* LS-DYNA (Finite Element (FE) program)

## Repository content
This repository contains following files.
* Bayesian_optimization_main_script.m : Main script to launch the optimization
* next_sample.m : 
* squared_exponential_n_dim.m : squared exponential covariance function definition
* surrogate_model.m : Gaussian process surrogate model
* lower_confidence_bound.m : lower confidence bound acquisition function
* get_reactfrcFEA.m : script to read the FE reaction force 
* get_residual.m : script to compute the residual value between FE and experimental data
* mean_square_error.m : script to compute mean squared error
* GetDynaReplace.m : script to modify the LS-DYNA input file with appropriate material constants
* Test31824.k : Base LS-DYNA input file
* Material_Matrix_Total.mat : Initial 30 FE simulations
* Anterior_Mean_101_20P.mat : MATLAB processed experimental mean data for the anterior region of meniscus 
* Anterior_SD_101_20P.mat  : MATLAB processed experimental standard deviation data for the anterior region of meniscus
* Middle_Mean_101_20P.mat : MATLAB processed experimental mean data for the middle region of meniscus
* Middle_SD_101_20P.mat : MATLAB processed experimental standard deviation data for the middle region of meniscus
* Posterior_Mean_101_20P.mat : MATLAB processed experimental mean data for the posterior region of meniscus
* Posterior_SD_101_20P.mat : MATLAB processed experimental standard deviation data for the posterior region of meniscus

## Usage
The provided results can be reproduced by following steps:
1. Download all the repository content to a directory
2. Give appropriate filepath to the LS-DYNA .exe in 'Bayesian_optimization_main_script.m' file
3. Launch 'Bayesian_optimization_main_script.m' using MATLAB


