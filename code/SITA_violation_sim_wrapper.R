#' ARGUMENTS:
#' rho (float) describes angle between propensity and prognostic model
#' p (int) number of covariates - must be at least 2
#' nsim (int) number of simulations

source("basic_sim_functions.R")
source("fullmatch_sim_functions.R")
source("pairmatch_sim_functions.R")
source("SITA_violation_sim_functions.R")

out_file <- "xSITA_results_"

# get user arguments
args <- commandArgs(trailingOnly=TRUE)

# parse args
rho <- as.numeric(args[1])
p <- as.numeric(args[2])
nsim <- as.numeric(args[3])

# defaults
true_mu <- "X1/3 - 3"
sigma <- 1
tau <- 1
ks <- c(1, 2, 3, 5)
N <- 2000
full <- F
nu <- 0.1
prop_model <- formula(t ~ . - mu - y - U) 
prog_model <- formula(y ~ . - mu - t - U)

run_sim <- function(rho = 0.1, p = 10, nsim = 10,
                    out_file = "test_", true_mu = "X1/3 - 3", nu = 0.1,
                    ks = 1:10, sigma = 1, tau = 1, N = 2000,
                    prog_model, prop_model,
                    full = FALSE) {
  t1 <- proc.time()
  
  message("********************")
  message("Simulation parameters:")
  message(paste("N:", N))
  message(paste("Rho:", rho))
  message(paste("p:",p))
  message(paste("nsim:", nsim))
  message(paste("true_mu:", true_mu))
  message(paste("nu:", nu))
  message(paste("sigma:", sigma))
  message(paste("tau:", tau))
  message(paste("full:", full))
  message("********************")
  
  # simulate
  if (full){
    exit("fullmatching with SITA violation is not yet supported")
    results <- replicate(nsim, simulate_fullmatch(generate_data(N=N, rho=rho,
                                                                p = p,
                                                                true_mu = true_mu,
                                                                sigma = sigma,
                                                                tau = tau),
                                      verbose = TRUE),
                       simplify = FALSE) %>% 
      bind_rows()
  } else {
    results <- replicate(nsim,
                         simulate_pairmatch(generate_xSITA_data(N=N, rho=rho, p = p,
                                                                true_mu = true_mu, nu = nu,
                                                                sigma = sigma, tau = tau),
                                            prop_model = prop_model,
                                            prog_model = progn_model,
                                            verbose = TRUE, ks = ks),
                         simplify = FALSE) %>% 
      bind_rows()
  }
  
  message("********************")
  message("Simulations complete:")
  message("********************")
  
  # write to file
  if (!is.null(out_file)){
    write.csv(results, file = paste0(out_file, rho*10, "_", p, "_", nsim, "_nu_", nu*10), row.names = FALSE)
    message(paste("output file:", paste0(out_file, rho*10, "_", p, "_", nsim, "_nu_", nu)))
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
        sigma = sigma, tau = tau, ks = ks, N = N,
        prog_model = prog_model, prop_model = prop_model,
        full = full)
