# Data

## File naming convention

A simulation of NSIM replications of data with p covariate dimensions with correllation rho between prognosis and treatment would be saved as:

`angle_sigma1_results_rho_p_NSIM`


## log files

Each simulation has an associated `.log` file in the same directory.  The naming convention for these files is:

`R_sim_YYYY-MM-DD.log`

where R is 10 times the rho parameter for the assiociated simulation, and YYYY-MM-DD gives the date the simulations were started.  Later log files contain all of the parameters for the simulation, a series of output lines that are printed whenever a simulation is completed, and (at the end of the file) the time elapsed while this simulation was run.  Earlier log files are not as detailed; the parameters may be incomplete, and will be printed at the end of the file rather than the beginning.

## How to produce these files

See the instructions in the `code` directory.

