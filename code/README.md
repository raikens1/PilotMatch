# Buffalo - code

This subdirectory contains all of the code to run the simulations analyzed in this manuscript.

## Quick Start

Here are some quick steps to get things up and running. The process below will let you run pairmatching or fullmatching simulations.

0. (Optional) start AWS EC2 instance `Rocky-Buffalo` or `Rocky-2Buffalo`

1. Check/modify the parameters in `basic_sim_wrapper.R`
   - For defaults, see the section below titled "Simulation Parameter Notes"

2. Run `sh basic_sim_batch_driver.sh NSIM`, where `NSIM` is the number of simulations you'd like to run 
   - usually `NSIM=1000` for pairmatching and `NSIM=2000` for fullmatching).
   
The commands above will produce: 11 log files, and 11 result files.  Each file is numbered 0,....10, and contains simulations for different values of the parameter rho (for rho = 0, 0.1,... 1.0).  The log file will show you the simulation parameters and document any errors; the results file should give the estimate and gamma sensitivity for Mahalanobis, Propensity, and Buffalo matching for each simulation.

Once the data are obtained, they should be moved to an appropriate location in the `data` subdirectory. Sample code to analyze new data can be found in the `analysis` directory.

**Note:** If you are trying to measure the distance between matched pairs in your simulation, the process is slightly different.  In step 1, you should modify the parameters in `distance_sim_wrapper.R`.  In step 2, you should run `sh distance_sim_batch_driver.sh`.  These simulations tend to be faster (~3 hours maybe) because they don't sweep out 11 different values of `k` and `rho`.  Also note that distance simulations do not yet support fullmatching.

## File specifics

In specific, the files to run most simulations are:
- `basic_sim_functions.R` - based on `simulation_functions.R` by Dylan Greaves.  This file contains the functions to do pairmatching for most simulations.  This is also useful for making visualizations like in Figure 1.
- `fullmatch_sim_functions.R` - additional functions for fullmatching simulations.  Has `basic_sim_functions.R` as a dependency
- `basic_sim_wrapper.R` - general wrapper used to run batches of NSIM simulations with specified parameters. Can be run from the command line but is more often called by `basic_sim_batch_driver.sh`.  Has `basic_sim_functions.R` and `fullmatch_sim_functions.R` as dependencies.
- `basic_sim_batch_driver.sh` - a very simple shell script that calls `basic_sim_wrapper.R` 11 times to run simulations across each value of rho (0.0, 0.1, ... 1.0).  Each call to `basic_sim_wrapper.R` runs in the background in parallel and produces a `.log` file and a results file (`angle_sigma1_results...`).

The files to run distance simulations are:
- `distance_sim_functions.R` - functions to do pairmatchings and record the distance between matched pairs by different metrics.  Has `basic_sim_functions.R` as a dependency.
- `distance_sim_wrapper.R` - similar to `basic_sim_wrapper.R` except that it runs the functions from `distance_sim_functions.R` and it has slightly different default parameters.
- `distance_sim_batch_driver.sh` - similar to `basic_sim_batch_driver.sh`. Instead of calling the wrapper to run simulations over 11 values of rho, this script calls the wrapper to run simulations over 6 different sample sizes (`n = 2000, 1800, ... 1000`)

In practice, I never actually used the wrapper and batch driver for the distance sims that much.  Most of what I did was smaller scale simulations with fixed sample size, just for the purpose of diagnosis (see the relevant file in `analyses`).


## Simulation Parameter Notes
This section details some simulation parameter defaults

### Default parameters

For most simulations, the basic parameters are:

```
true_mu = "X1/3 - 3"
sigma = 1
tau = 1
ks = 1:10 (if applicable)
N = 2000
NSIM = 1000 for pairmatching (full = F) and 2000 for fullmatching (full = T)
rho = 0, 0.1, ... 1
p = 10
```

**Notes:** 
- It's worth noting that `mu` is actually not the notation we use for this parameter in the paper.  It's phi.
- For distance simulations, you will only use one rho and one k.  I have generally used `rho = 0.5` and `k = 3`.  This isn't set in stone, however, since none of these simulations made it into the final manuscript anyway.
- The number of simulations for fullmatch runs was increased because I found that the fullmatch results were especially noisy.  Even at 2000 simulations, estimates of the standard deviation can be noisy.

### Varying sample size

In order to decrease the sample size while keeping the number of treated individuals at about 100, I have had to tinker with the constant in the parameter mu.

```
 N       mu
 2000    "X1/3 - 3"
 1800    "X1/3 - 2.9"
 1600    "X1/3 - 2.75"
 1400    "X1/3 - 2.6"
 1200    "X1/3 - 2.45"
 1000    "X1/3 - 2.25"
```

### Quick and Dirty Simulations

I have run some simulations where I just wanted to get a quick idea of what would happen.  For this, I let `ks = c(1,4,8)`, and set `NSIM` to be something smaller (say, 100).  Simulations like these take about 15 minutes, but the results are much more noisy.
