#' ARGUEMENTS:
#' rho (float) describes angle between propensity and prognostic model
#' p (int) number of covariates - must be at least 2
#' nsim (int) number of simulations
library(pbapply)
source("../code/simulation_functions.R")

# set default arguements
rho <- 0.1
p <- 10
nsim <- 100
out_file <- "angle_sigma1_results_"


# get user arguements
args <- commandArgs(trailingOnly=TRUE)

# parse args
rho <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim <- as.numeric(args[3])

# simulate and write to file
results <- pbreplicate(nsim, simulate(generate_data(rho=rho, p = p)),
                     simplify = FALSE) %>% 
  bind_rows
write.csv(results, file = paste0(out_file, rho*10, "_", p, "_", nsim), row.names = FALSE)

message("********************")
message("Simulation complete:")
message(paste("Rho:", rho))
message(paste("p:",p))
message(paste("nsim:", nsim))
message(paste("output file:", out_file))
message("********************")pbapplu