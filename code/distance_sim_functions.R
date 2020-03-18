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

source("../code/basic_sim_functions.R")

#' @title Simulate for distances
#' @description perform matchings like simulate, but with fixed k, and return distances between matches
#'  rather than effect estimate and gamma 
#' @param df, a data.frame from generate_data
#' @param prop_model, the propensity score model
#' @param prog_model, the prognostic score model
#' @param k, the number of controls to match to each treated
#' @param true_rho the actual value of rho
#' @param gamma boolean for whether or not to compute gamma (takes longer if TRUE)
#' @return a data.frame with results from propensity, mahalanobis, and buffalo matching
simulate_for_distances <- function(df, 
                                   prop_model = formula(t ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10), 
                                   prog_model = formula(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10),
                                   verbose = FALSE,
                                   k = 3,
                                   true_rho = 0,
                                   gamma = FALSE){
  # if not enough controls, do less than 1:k, but print an error
  if(sum(df$t)*5 > nrow(df)){
    kmax = floor(nrow(df)/sum(df$t))-2
    k = min(k, kmax)
    message(paste0("Insufficient controls.  Doing 1:", k, " matching instead"))
  }
  
  # propensity score matching for k = 1:10
  propensity <- glm(prop_model, family = binomial(), data = df)
  
  prop_match <- pairmatch(propensity, controls = k, df)
  
  # 1:2 mahalanobis matching to select data to use for prognostic model
  mahal_dist <- match_on(prop_model, method = "mahalanobis", data = df)
  mahal_match <- pairmatch(mahal_dist, controls = 2, df) 
  
  # perform prognostic score matching for k
  prog_list <- prognostic_match(df, propensity, mahal_match, prog_model, k, params_only = TRUE)
  
  # mahalanobis matching for k
  m_match <- pairmatch(mahal_dist, controls = k, df)
  
  prop_df <- data_frame(method = "propensity",
                        k = k,
                        mean_dist = true_mean_dist(df, prop_match, true_rho),
                        var_dist = true_var_dist(df, prop_match, true_rho),
                        emp_mean_dist = emp_mean_dist(df, prop_match, propensity, prog_list$prognostic),
                        prog_dist = prog_dist(df, prop_match, true_rho),
                        prop_dist = prop_dist(df, prop_match),
                        estimate = att_estimate(reformat(df, prop_match, k)),
                        gamma = ifelse(gamma, gamma_sensitivity(reformat(df, prop_match, k)), NA))
  
  prog_df <- data_frame(method = "prognostic",
                        k = k,
                        mean_dist = true_mean_dist(prog_list$df, prog_list$prog_match, true_rho),
                        var_dist = true_var_dist(prog_list$df, prog_list$prog_match, true_rho),
                        emp_mean_dist = emp_mean_dist(prog_list$df, prog_list$prog_match, propensity, prog_list$prognostic),
                        prog_dist = prog_dist(prog_list$df, prog_list$prog_match, true_rho),
                        prop_dist = prop_dist(prog_list$df, prog_list$prog_match),
                        estimate = att_estimate(reformat(prog_list$df, prog_list$prog_match, k)),
                        gamma = ifelse(gamma, gamma_sensitivity(reformat(prog_list$df, prog_list$prog_match, k)), NA))
  
  mahal_df <- data_frame(method = "mahalanobis",
                         k = k,
                         mean_dist = true_mean_dist(df, m_match, true_rho),
                         var_dist = true_var_dist(df, m_match, true_rho),
                         emp_mean_dist = emp_mean_dist(df, m_match, propensity, prog_list$prognostic),
                         prog_dist = prog_dist(df, m_match, true_rho),
                         prop_dist = prop_dist(df, m_match),
                         estimate = att_estimate(reformat(df, m_match, k)),
                         gamma = ifelse(gamma, gamma_sensitivity(reformat(df, m_match, k)), NA))
  if (verbose){
    message("Completed One Simulation")
  }
  # return results for prop, prog, and mahal
  return(bind_rows(prop_df, prog_df, mahal_df))
}

#' @title Get True Distances
#' @description Return the mean mahalanobis distance between matched t and c individuals in terms of true scores
#' @param df data.frame of all individuals
#' @param match the result of a call to fullmatch or pairmatch
#' @param true_rho (float) the true value of rho
#' @return (float) mean mahalanobis distance between matched t and c individuals from match
true_mean_dist <- function(df, match, true_rho){
  ndf <- df %>% mutate( psi = true_rho*X1 + sqrt(1-true_rho^2)*X2) 
  dists <- match_on(t ~ psi + mu, data = ndf)
  dist_list <- matched.distances(match, dists)
  return(mean(unlist(dist_list)))
}

#' @title Get variance in True Distances
#' @description Return the mean mahalanobis distance between matched t and c individuals in terms of true scores
#' @param df data.frame of all individuals
#' @param match the result of a call to fullmatch or pairmatch
#' @param true_rho (float) the true value of rho
#' @return (float) variance in mahalanobis distance between matched t and c individuals from match
true_var_dist <- function(df, match, true_rho){
  df <- df %>% mutate( psi = true_rho*X1 + sqrt(1-true_rho^2)*X2) 
  dists <- match_on(t ~ psi + mu, data = df)
  dist_list <- matched.distances(match, dists)
  return(var(unlist(dist_list)))
}

#' @title Get Empirical Distances
#' @description Return the mean mahalanobis distance between matched t and c individuals in terms of estimated scores
#' @param df data.frame of all individuals
#' @param match the result of a call to fullmatch or pairmatch
#' @param true_rho (float) the true value of rho
#' @return (float) mean mahalanobis distance between matched t and c individuals from match
emp_mean_dist <- function(df, match, propensity, prognosis){
  df <- df %>% mutate(prog = predict(prognosis, df), prop = predict(propensity, df)) 
  dists <- match_on(t ~ prog + prop, data = df)
  dist_list <- matched.distances(match, dists)
  return(mean(unlist(dist_list)))
}

prog_dist <- function(df, match, true_rho){
  ndf <- df %>% mutate( psi = true_rho*X1 + sqrt(1-true_rho^2)*X2, m = as.character(match)) 
  treat_df <- ndf %>% filter(t == 1, !is.na(m))  %>% mutate(psi_treat = psi) %>% select(c(psi_treat, m))
  control_df <- ndf %>% filter(t == 0, !is.na(m)) %>% full_join(treat_df, by = "m") %>% mutate(prog_dist = abs(psi-psi_treat))
  return(mean(control_df$prog_dist))
}

prop_dist <- function(df, match, true_rho){
  ndf <- df %>% mutate(m = as.character(match)) 
  treat_df <- ndf %>% filter(t == 1, !is.na(m))  %>% mutate(mu_treat = mu) %>% select(c(mu_treat, m))
  control_df <- ndf %>% filter(t == 0, !is.na(m)) %>% full_join(treat_df, by = "m") %>% mutate(prop_dist = abs(mu-mu_treat))
  return(mean(control_df$prop_dist))
}