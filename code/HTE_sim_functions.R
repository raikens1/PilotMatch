require("optmatch", quietly = T)
require(dplyr, quietly = T)
require(magrittr, quietly = T)
require(ggplot2, quietly = T)
require(rlang, quietly = T)
require(tidyr, quietly = T)
require(sensitivitymw, quietly = T)
require(tidyselect, quietly = T)
require(bindrcpp, quietly = T)
require(sensitivityfull, quietly = T)

#' @title Generate simulated data with treatment effect heterogeneity
#' @description Data set contains measured covariates X, outcome Y, treatment
#'   assignment t, treatment effect tau, and logit(propensity), mu.
#'
#'   covariate data ~ normal(0,1)
#'   mu = true_mu
#'   t ~ binom(p = 1/(1+exp(-mu)))
#'   y ~ rho * X1 + sqrt(1-rho^2) * X2 + tau * t + epsilon
#'   epsilon ~ normal(0, 1)
#'
#' @param N numeric, sample size
#' @param p numeric, number of features
#' @param true_mu string formula giving true propensity score linear model
#' @param rho numeric between 0 and 1.  0 => prog orthogonal to prop, 1=> prog
#'   || prop
#' @param sigma numeric noise to be added to y. y += sigma*rnorm(0,1)
#' @param true_tau treatment effect
#' @return data.frame of covariates, y, t, and mu
generate_data_HTE <- function(N = 2000,
                          p = 10,
                          true_mu = "X1/3-3", 
                          rho = 0,
                          sigma = 1,
                          true_tau = "1 + rnorm(N, sd = 0.25)") {
  df <- data.frame(matrix(rnorm(p*N), ncol = p))
  df <- df %>% mutate(mu = !!parse_quosure(true_mu))
  df <- df %>% mutate(tau = !!parse_quosure(true_tau))
  df <- df %>% mutate(t = rbinom(n = N, size = 1, prob = 1/(1+exp(-mu))))
  df <- df %>% mutate(y = tau*t + rho*X1 + sqrt(1-rho^2)*X2)
  noise <- rnorm(N)
  df$y <- df$y + sigma*noise
  return(df)
}

get_SATT <- function(df){
  treated_df <- df %>% filter(t == 1)
  
  return(mean(treated_df$tau))
}