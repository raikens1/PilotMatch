# Simulation Parameter Notes
This document details some simulation parameter defaults

## Default parameters

For most simulations, the basic parameters are:
true_mu = "X1/3 - 3"
sigma = 1
tau = 1
ks = 1:10 (if applicable)
N = 2000
NSIM = 1000 for pairmatching and 2000 for fullmatching
rho = 0, 0.1, ... 1
p = 10

## Varying sample size

When N = 2000, true_mu = "X1/3 - 3" (default)
When N = 1600, true_mu = "X1/3 - 2.75"
When N = 1100, true_mu = "X1/3 - 2.35"
When N = 1000, true_mu = "X1/3 - 2.25"
