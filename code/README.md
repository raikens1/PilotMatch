# Buffalo - code

This subdirectory contains all of the code to run the simulations analyzed in this manuscript.

## Quick Start

Here are some quick steps to get things up and running. The process below will let you run pairmatching or fullmatching simulations.

0. (Optional) start AWS EC2 instance `Rocky-Buffalo` or `Rocky-2Buffalo`

1. Check/modify the parameters in basic_sim_wrapper.sh_ 
   - For defaults, see the section below titled "Simulation Parameter Notes"

2. Run `sh basic_sim_batch_driver.sh NSIM`, where  



## Simulation Parameter Notes
This section details some simulation parameter defaults

### Default parameters

For most simulations, the basic parameters are:
true_mu = "X1/3 - 3"
sigma = 1
tau = 1
ks = 1:10 (if applicable)
N = 2000
NSIM = 1000 for pairmatching and 2000 for fullmatching
rho = 0, 0.1, ... 1
p = 10

### Varying sample size

 N       mu
 2000    "X1/3 - 3"
 1800    "X1/3 - 2.9"
 1600    "X1/3 - 2.75"
 1400    "X1/3 - 2.6"
 1200    "X1/3 - 2.45"
 1000    "X1/3 - 2.25"
