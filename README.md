# Buffalo
*Authors:* Rocky Aikens, Dylan Greaves, Michael Baiocchi

This repository contains the code and data used to write an upcoming paper on using the prognostic score in matching for causal inference from observational studies

## Contents of this repo

- **Analyses** - all vizualizations, interpretations, and analysis of simulated data are in .Rmd files here. The most important files are `Main_Figures.Rmd`, which makes all of the main figures for the manuscript and `Supplementary_Figures.Rmd`, which contains all the supplementary figures.
- **Code** - all code (R and bash) to run the simulations is here.
- **Data** - all data produced from simulations is here.
- **Figures** - current versions of all figures are here.  These figures are automatically produced from the code in `analyses/Main_Figures.Rmd`

**Code** and **Data** each contain a README further explaining the contents of that subdirectory.

## Manuscript

The manuscript "Using the Prognostic Score to Reduce Heterogeneity in Observational Studies" (Aikens, Greaves, and Baiocchi) is written but not yet submitted.  It is written in a shared overleaf document owned by MB.

## Other note: AWS

Much of the computational heavylifting to produce these simulations was done on an Amazon Web Server EC2 instance.  There are two EC2 instances which I have made for this project: `Rocky-Buffalo` and `Rocky-2Buffalo`.  On each of these instances, a standard batch of 1000 simulations can take 6-24 hours, depending on the sample size and specific simulation parameters.  Pairmatching tends to take longer, Fullmatching tends to be faster.



