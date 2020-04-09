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
require(glmnet, quietly = T)

source("../code/basic_sim_functions.R")

#' @title simulate Pairmatch with lasso fit propensity and prognostic models
#' @param df, a data.frame from generate_data
#' @param prop_model, the propensity score model
#' @param prog_model, the prognostic score model
#' @param mahal_model, the mahalanobis distance matching model (t ~ covariates
#'   to match on)
#' @param ks, a vector of positive numbers indicating the number of controls to
#'   match to each treated
#' @return a data.frame with results from propensity, mahalanobis, and buffalo
#'   matching
simulate_pairmatch_lasso <- function(df, 
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
  
  # build lasso propensity score model
  x_all <- df %>%
    select_(.dots = labels(terms(prop_model))) %>%
    as.matrix()
  
  y_all <- df %>%
    select_(.dots = all.vars(prop_model)[1]) %>%
    as.matrix()
  
  cvprop <- cv.glmnet(x_all, y_all, family = "binomial")
  prop_scores <- predict(cvprop, newx = x_all, s = "lambda.min")
  
  f <- function(k){ 
    prop_match <- pairmatch(t ~ propscore,
                            controls = k,
                            data = mutate(df, propscore = prop_scores))
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
    prognostic_match_lasso(df = df, prop_scores = prop_scores, 
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


#' @title prognostic_match_lasso
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
prognostic_match_lasso <- function(df, prop_scores, match_assignment,
                             prog_model, k, params_only = FALSE) {
  
  df$m <- match_assignment
  df$row <- 1:nrow(df)
  df$propscore <- as.numeric(prop_scores)
  
  # select pilot set
  pilot_set <- df %>% 
    filter(!is.na(m)) %>%
    filter(t==0) %>%
    group_by(m) %>%
    sample_n(size = 1) %>%
    ungroup()
  
  # fit lasso prognostic model
  x_pilot <- pilot_set %>%
    select_(.dots = labels(terms(prog_model))) %>%
    as.matrix()
  
  y_pilot <- pilot_set %>%
    select_(.dots = all.vars(prog_model)[1]) %>%
    as.matrix()
  
  cvprog <- cv.glmnet(x_pilot, y_pilot)
  
  # match analysis set
  analysis_set <- df[-pilot_set$row, ]
  
  x_analysis <- analysis_set %>%
    select_(.dots = labels(terms(prog_model))) %>%
    as.matrix()
  
  prog_scores <- predict(cvprog, newx = x_analysis, s = "lambda.min")
  
  analysis_set <- analysis_set %>% 
    mutate(progscore = prog_scores) 
  prog_dist <- match_on(t ~ progscore + propscore, data = analysis_set)
  prog_match <- pairmatch(prog_dist, controls = k, data = analysis_set) 
  
  # return
  if(params_only){
    return(list(df = analysis_set, prog_match = prog_match, prognostic = prognostic))
  } else{
    return(reformat(analysis_set, prog_match, k))}
}
