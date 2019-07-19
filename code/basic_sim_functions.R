require("optmatch")
require(dplyr)
require(magrittr)
require(ggplot2)
require(rlang)
require(tidyr)
require(sensitivitymw)
require(tidyselect)
require(bindrcpp)
require(sensitivityfull)

#' @title Generate simulated data: X, y, t, and mu
#' @description covariate data ~ normal(0,1); mu = true_mu; 
#'              t ~ binom(p = 1/(1+exp(-mu))); y ~ X1 + X2 + t
#' @param N numeric, sample size
#' @param p numeric, number of features
#' @param true_mu string formula giving true propensity score linear model
#' @param rho numeric between 0 and 1.  0 => prog orthogonal to prop, 1=> prog || prop
#' @param sigma numeric noise to be added to y. y += sigma*rnorm(0,1)
#' @param tau numeric additive treatment effect
#' @return data.frame of covariates, y, t, and mu
generate_data <- function(N = 2000,
                          p = 10,
                          true_mu = "X1/3-3", 
            			        rho = 0,
                          sigma = 1,
            			        tau = 1) {
  df <- data.frame(matrix(rnorm(p*N), ncol = p))
  df <- df %>% mutate(mu = !!parse_quosure(true_mu))
  df <- df %>% mutate(t = rbinom(n = N, size = 1, prob = 1/(1+exp(-mu))))
  df <- df %>% mutate(y = tau*t + rho*X1 + sqrt(1-rho^2)*X2)
  noise <- rnorm(N)
  df$y <- df$y + sigma*noise
  return(df)
}

#' @title Reformat Data based on a matching
#' @param df an R data.frame of data from generate_data()
#' @param match optmatch object giving a matching
#' @param n_control int - number of controls matched to each treated
#' @return a new data.frame of y values where each row is a matched set
#'  and columns are individuals from the set (1st col treated, all following cols control)
reformat <- function(df, match, n_control) {
  ndf <- df %>% 
    mutate(m = match) %>%
    select(y,t,m) %>%
    filter(!is.na(m)) 
  
  ndf <- ndf %>% 
    arrange(m, desc(t)) %>% 
    mutate(id = rep(1:(n_control + 1), sum(ndf$t))) %>%
    select(-t) %>%
    spread(id, y) %>% 
    select(-m)

  return(ndf)
}

#' @title  att_estimate
#' @param ndf a data.frame of y values by matched sets
#' @returns att, numeric
att_estimate <- function(ndf) {
  return(senmwCI(ndf, method = "t")$PointEstimate[[1]])
}

#' @title gamma_sensitivity
#' @param ndf a data.frame of y values by matched sets
#' @returns gamma sensitivity for data
gamma_sensitivity <- function(ndf) {
  delta <- 0.01
  g <- 1 + delta
  while (g < 200){
    pval <- senmw(ndf, gamma = g, method = "t")$pval
    if (pval > 0.05){
	return(g - delta)
    }
    g <- g + delta
  }
}

#' @title prognostic_match
#' @param df a data.frame from generate_data()
#' @param propensity glm object for fitted propensity score model
#' @param match_assignment optmatch object with assignment from mahalanobis 1:2 matching
#' @param prog_model formula for prognostic model
#' @param n_control number of control individuals to match to each treated
#' @param params_only (boolean) if TRUE, return not_selected and match object, rather than reformatted data frame
#' @return a data.frame of y reformatted by matching assignment according to buffalo method
prognostic_match <- function(df, propensity, match_assignment, prog_model, n_control, params_only = FALSE) {
  df$m <- match_assignment
  df$row <- 1:nrow(df)
  n_t<- sum(df$t)

  selected <- df %>% 
    filter(!is.na(m)) %>%
    filter(t==0) %>%
    group_by(m) %>%
    sample_n(size = 1) %>%
    ungroup()
  
  prognostic <- lm(y ~ . - mu - t, data = dplyr::select(selected, -c(row, m)))
  not_selected <- df[-selected$row, ]
  not_selected <- not_selected %>% 
			mutate(progscore = predict(prognostic, not_selected)) %>%
			mutate(propscore = predict(propensity, not_selected))
  prog_dist <- match_on(t ~ progscore + propscore, data = not_selected)
  prog_match <- pairmatch(prog_dist, controls = n_control, data = not_selected) 
  
  if(params_only){
    return(list(df = not_selected, prog_match = prog_match, prognostic = prognostic))
  } else{
  return(reformat(not_selected, prog_match, n_control))}
}