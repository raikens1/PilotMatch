#' ARGUEMENTS:
#' rho (float) describes angle between propensity and prognostic model
#' p (int) number of covariates - must be at least 2
#' nsim (int) number of simulations

source("simulation_functions.R")

out_file <- "dist_results_"

# Parameters N and mu for different N
# N       mu
# 2000    "X1/3 - 3"
# 1800    "X1/3 - 2.9"
# 1600    "X1/3 - 2.75"
# 1400    "X1/3 - 2.6"
# 1200    "X1/3 - 2.45"
# 1000    "X1/3 - 2.25"

Ns <- c(2000, 1800, 1600, 1400, 1200, 1000)
mus <- c("X1/3 - 3", "X1/3 - 2.9", "X1/3 - 2.75", "X1/3 - 2.6", "X1/3 - 2.45", "X1/3 - 2.25")

# defaults
true_mu <- "X1/3 - 3"
sigma <- 1
tau <- 1
k <- 3
N <- 2000
full <- F
p <- 10

# get user arguements
args <- commandArgs(trailingOnly=TRUE)

# parse args
rho <- as.numeric(args[1])
N <- as.numeric(args[2])
nsim <- as.numeric(args[3])
true_mu <- mus[which(Ns == N)]

run_sim <- function(rho = 0.1, p = 10, nsim = 10,
                    out_file = "test_", true_mu = "X1/3 - 3", 
                    k = 3, sigma = 1, tau = 1, N = 2000,
                    full = FALSE) {
  t1 <- proc.time()
  
  message("********************")
  message("Simulation parameters:")
  message(paste("N:", N))
  message(paste("Rho:", rho))
  message(paste("p:",p))
  message(paste("nsim:", nsim))
  message(paste("true_mu:", true_mu))
  message(paste("sigma:", sigma))
  message(paste("tau:", tau))
  message(paste("full:", full))
  message(paste("k:", k))
  message("********************")
  
  # simulate
  if (full){
    message("Full matching not yet supported. Exiting.")
  } else {
    results <- replicate(nsim, simulate_for_distances(generate_data(N=N, rho=rho, p = p, true_mu = true_mu, sigma = sigma, tau = tau),
                                        verbose = TRUE, k = k, true_rho = rho),
                         simplify = FALSE) %>% 
      bind_rows()
  }
  
  message("********************")
  message("Simulations complete:")
  message("********************")
  
  # write to file
  if (!is.null(out_file)){
    write.csv(results, file = paste0(out_file, rho*10, "_", N, "_", nsim), row.names = FALSE)
    message(paste("output file:", paste0(out_file, rho*10, "_", N, "_", nsim)))
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

run_sim(rho = rho, p = p, nsim = nsim, 
        out_file = out_file, true_mu = true_mu, 
        sigma = sigma, tau = tau, k = k, N = N,
        full = full)

print(true_mu)