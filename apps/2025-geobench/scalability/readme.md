### Scalability studies

This folder contains a simple app to solve the Stokes system with the FGMRES solver with multigrid preconditioning. In default setting the runtimes are reported for 5 FGMRES iterations. Other solver paramters such as number of smoothing steps, coarse grid solver iterations, etc can be adjusted in the corresponding parameter `.prm` file. The configuration of the setup can also be adjusted in the `.prm` file. The `generate_prm_runner_scaling.py` file can be used to generate folders, runner and parameter files to run multiple scalability studies.

