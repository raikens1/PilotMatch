# Data

## File naming convention

A simulation of n replications of data with p covariate dimensions as part of experiment q would be documented as:

"angle_sigma1_results_n_p_q"

## How to produce these files

A basic command to produce files like this might be:

for i in {1..50}; do Rscript ../R/prognostic_angle.R ${i} 10; done

The first argument after the script name being the experiment number (q above), and the second being the covariate dimension (p above)

