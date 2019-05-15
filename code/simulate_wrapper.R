#' ARGUEMENTS:
#' rho (float) describes angle between propensity and prognostic model
#' p (int) number of covariates - must be at least 2
#' nsim (int) number of simulations


source("simulation_functions.R")

out_file <- "angle_sigma1_results_"

# get user arguements
args <- commandArgs(trailingOnly=TRUE)

# parse args
rho <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim <- as.numeric(args[3])

# defaults
true_mu <- "X1 - 10/3"
sigma <- 1
tau <- 0.5
ks <- c(1,4,8)

run_sim <- function(rho = 0.1, p = 10, nsim = 10, out_file = "test_", true_mu = true_mu, ks = 1:10, sigma = sigma, tau = tau, N = 1000) {
  t1 <- proc.time()
  # simulate
  results <- replicate(nsim, simulate(generate_data(N=N, rho=rho, p = p, true_mu = true_mu, sigma = sigma, tau = tau),
                                      verbose = TRUE, ks = ks),
                       simplify = FALSE) %>% 
    bind_rows()
  
  message("********************")
  message("Simulation complete:")
  message(paste("Rho:", rho))
  message(paste("p:",p))
  message(paste("nsim:", nsim))
  message(paste("true_mu", true_mu))
  message(paste("sigma", sigma))
  message(paste("tau", tau))
  
  # write to file
  if (!is.null(out_file)){
    write.csv(results, file = paste0(out_file, rho*10, "_", p, "_", nsim), row.names = FALSE)
    message(paste("output file:", paste0(out_file, rho*10, "_", p, "_", nsim)))
  } else {
    message("Results not saved to file.")
  }
  message("Time elapsed:")
  message(print(proc.time() - t1))
  message("********************")
  
  if (is.null(out_file)){
    return(results)
    } else{
      return(1)
    }

}

run_sim(rho = rho, p = p, nsim = nsim, out_file = out_file, true_mu = true_mu, sigma = sigma, tau = tau, ks = ks, N = 1000)
 