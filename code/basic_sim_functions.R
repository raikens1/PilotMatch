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

#' @title Generate simulated data with a SITA violation
#' @description Data set contains measured covariates X, outcome Y, treatment
#'   assignment t, and logit(propensity), mu.
#'
#'   covariate data ~ normal(0,1); mu = true_mu; t ~ binom(p = 1 / (1 +
#'   exp(-mu))); y ~ rho * X1 + sqrt(1-rho^2) * X2 + tau * t + epsilon epsilon ~
#'   normal(0, 1)
#'
#' @param N numeric, sample size
#' @param p numeric, number of features
#' @param true_mu string formula giving true propensity score linear model
#' @param rho numeric between 0 and 1.  0 => prog orthogonal to prop, 1=> prog
#'   || prop
#' @param sigma numeric noise to be added to y. y += sigma*rnorm(0,1)
#' @param tau numeric additive treatment effect
#' @return data.frame of covariates, y, t, and mu
generate_data <- function(N = 2000,
                          p = 10,
                          true_mu = "X1/3-3", 
            			        rho = 0,
                          sigma = 1,
            			        tau = 1) {
  
  df <- data.frame(matrix(rnorm(p*N), ncol = p)) %>%
    mutate(mu = !!parse_quosure(true_mu),
           t = rbinom(n = N, size = 1, prob = 1 / (1 + exp(-mu))),
           y = tau * t + rho * X1 + sqrt(1 - rho ^ 2) * X2 + rnorm(N, sd = sigma))

  return(df)
}

#' @title Reformat Data based on a matching
#' @param df an R data.frame of data from generate_data()
#' @param match optmatch object giving a matching
#' @param k int - number of controls matched to each treated
#' @return a new data.frame of y values where each row is a matched set
#'  and columns are individuals from the set (1st col treated, all following cols control)
reformat <- function(df, match, k) {
  ndf <- df %>% 
    mutate(m = match) %>%
    select(y,t,m) %>%
    filter(!is.na(m)) 
  
  ndf <- ndf %>% 
    arrange(m, desc(t)) %>% 
    mutate(id = rep(1:(k + 1), sum(ndf$t))) %>%
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
#' @param match_assignment optmatch object with assignment from mahalanobis 1:2
#'   matching
#' @param prog_model formula for prognostic model
#' @param k number of control individuals to match to each treated
#' @param params_only (boolean) if TRUE, return analysis_set and match object,
#'   rather than reformatted data frame
#' @return a data.frame of y reformatted by matching assignment according to
#'   buffalo method
prognostic_match <- function(df, propensity, match_assignment,
                             prog_model, k, params_only = FALSE) {
  
  df$m <- match_assignment
  df$row <- 1:nrow(df)

  # select pilot set
  pilot_set <- df %>% 
    filter(!is.na(m)) %>%
    filter(t==0) %>%
    group_by(m) %>%
    sample_n(size = 1) %>%
    ungroup()
  
  # fit prognostic model
  prognostic <- lm(prog_model, data = dplyr::select(pilot_set, -c(row, m)))
  
  # match analysis set
  analysis_set <- df[-pilot_set$row, ]
  
  analysis_set <- analysis_set %>% 
			mutate(progscore = predict(prognostic, analysis_set)) %>%
			mutate(propscore = predict(propensity, analysis_set))
  prog_dist <- match_on(t ~ progscore + propscore, data = analysis_set)
  prog_match <- pairmatch(prog_dist, controls = k, data = analysis_set) 
  
  # return
  if(params_only){
    return(list(df = analysis_set, prog_match = prog_match, prognostic = prognostic))
  } else{
  return(reformat(analysis_set, prog_match, k))}
}