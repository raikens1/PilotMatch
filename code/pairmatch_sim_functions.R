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

#' @title simulate Pairmatch
#' @param df, a data.frame from generate_data
#' @param prop_model, the propensity score model
#' @param prog_model, the prognostic score model
#' @param mahal_model, the mahalanobis distance matching model (t ~ covariates
#'   to match on)
#' @param ks, a vector of positive numbers indicating the number of controls to
#'   match to each treated
#' @return a data.frame with results from propensity, mahalanobis, and buffalo
#'   matching ??
simulate_pairmatch <- function(df, 
                     prop_model = formula(t ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10), 
                     prog_model = formula(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10),
                     mahal_model = formula(t ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10),
                     verbose = FALSE,
                     ks = 1:5, HTE = F) {
  
  # if not enough controls, do less than 1:k, but print an error
  if(sum(df$t)*12 > nrow(df)){
    kmax = floor(nrow(df)/sum(df$t))-2
    ks = 1:min(10, kmax)
    message(paste0("Insufficient controls.  Doing 1:1 to 1:", min(10, kmax), " matching instead"))
  }
  
  # propensity score matching 
  # store sensitivity and att in prop_df
  propensity <- glm(prop_model, family = binomial(), data = df)
  
  f <- function(k){ 
    prop_match <- pairmatch(propensity, controls = k, data = df)
    return(reformat(df, prop_match, k))
  }
  ndfs <- sapply(ks, f, simplify = FALSE)
  prop_df <- data_frame(method = "propensity",
                        k = ks,
                        estimate = sapply(ndfs, att_estimate),
                        gamma = sapply(ndfs, gamma_sensitivity))
  
  # 1:2 mahalanobis matching to select data to use for prognostic model
  mahal_dist <- match_on(mahal_model,
                         method = "mahalanobis", data = df)
  mahal_match <- pairmatch(mahal_dist, controls = 2, df) 
  
  # perform prognostic score matching for k in ks
  # store sensitivity and att in prog_df
  g <- function(k){
    prognostic_match(df = df, propensity = propensity, 
                     match_assignment = mahal_match,
                     prog_model =  prog_model, k = k)
  }
  
  ndfs <- sapply(ks, g, simplify = FALSE)
  prog_df <- data_frame(method = "prognostic",
                        k = ks,
                        estimate = sapply(ndfs, att_estimate),
                        gamma = sapply(ndfs, gamma_sensitivity))
  
  # mahalanobis matching for k = 1:10
  # store sensitivity and att in mahal_df
  h <- function(k){
    m_match <- pairmatch(mahal_dist, controls = k, df)
    return(reformat(df, m_match, k))
  }
  
  ndfs <- sapply(ks, h, simplify = FALSE)
  mahal_df <- data_frame(method = "mahalanobis",
                         k = ks,
                         estimate = sapply(ndfs, att_estimate),
                         gamma = sapply(ndfs, gamma_sensitivity))
  if (verbose){
    message("Completed One Simulation")
  }
  
  result <- bind_rows(prop_df, prog_df, mahal_df) 
  
  if (HTE == T){
    result$SATT = get_SATT(df)
  }
  
  # return results for prop, prog, and mahal
  return(result)
}
