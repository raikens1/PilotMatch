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

#' @title simulate full match
#' @param df, a data.frame from generate_data
#' @param prop_model, the propensity score model
#' @param prog_model, the prognostic score model
#' @param ks, a vector of positive numbers indicating the number of controls to match to each treated
#' @return a data.frame with results from propensity, mahalanobis, and buffalo matching
simulate_fullmatch <- function(df, 
                               prop_model = formula(t ~ . - mu - y), 
                               prog_model = formula(y ~ . - mu - t),
                               verbose = FALSE) {
  
  # propensity score full matching
  propensity <- glm(prop_model, family = binomial(), data = df)
  
  prop_match <- fullmatch(propensity, data = df)
  prop_data <- df %>% mutate(subclass = as.character(prop_match))
  prop_df <- data_frame(method = "propensity",
                        estimate = fullmatch_estimate(prop_data))
  
  # 1:2 mahalanobis matching to select data to use for prognostic model
  mahal_dist <- match_on(prop_model, method = "mahalanobis", data = df)
  mahal_match <- pairmatch(mahal_dist, controls = 2, df) 
  
  # perform prognostic score full matching
  prog_data <- prognostic_fullmatch(df, propensity, mahal_match, prog_model)
  prog_df <- data_frame(method = "prognostic",
                        estimate = fullmatch_estimate(prog_data))
  
  # perform mahalanobis full matching
  m_match <- fullmatch(mahal_dist, data = df)
  m_data <- df %>% mutate(subclass = as.character(m_match))
  mahal_df <- data_frame(method = "mahalanobis",
                         estimate = fullmatch_estimate(m_data))
  if (verbose){
    message("Completed One Simulation")
  }
  # return results for prop, prog, and mahal
  return(bind_rows(prop_df, prog_df, mahal_df))
}


#' @title prognostic full match
#' @param df a data.frame from generate_data()
#' @param propensity glm object for fitted propensity score model
#' @param match_assignment optmatch object with assignment from mahalanobis 1:2 matching
#' @param prog_model formula for prognostic model
#' @param n_control number of control individuals to match to each treated
#' @param params_only (boolean) if TRUE, return not_selected and match object, rather than reformatted data frame
#' @return a data.frame of y reformatted by matching assignment according to buffalo method
prognostic_fullmatch <- function(df, propensity, match_assignment, prog_model, params_only = FALSE) {
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
  prog_match <- fullmatch(prog_dist, data = not_selected) 
  
  
  if(params_only){
    return(list(df = not_selected, prog_match = prog_match))
  } else{
    return(not_selected %>% mutate(subclass = as.character(prog_match)))}
}

#' @title Fullmatch effect estimate
#' @description returns effect estimate from fullmatch using sensitivityfull
#' @param match_data data.frame with subclass column giving fullmatch assignments
#' @return float effect estimate from fullmatch
fullmatch_estimate <- function(match_data){
  J <- max(table(match_data$subclass))
  I <- length(unique(match_data$subclass))
  
  outcome_df <- match_data %>% group_by(subclass) %>% do(fullmatch_reformat(., J))
  
  y <- matrix(outcome_df$y, I, J, byrow = TRUE)
  
  treated <- !(outcome_df %>% group_by(subclass) %>% summarize(flip = first(flip)))$flip
  
  return(senfmCI(y = y, treated1 = treated)$PointEstimates[1])
}

#' @title Fullmatch Reformat
#' @description helper function to reformat the data from a fullmatch; called by fullmatch_estimate
#' @param sub_df data.frame with the individuals in a given subclass
#' @param J int the size of the largest matched set
#' @return new_df, a reformatted J by 4 data frame with rows of NAs for filler
fullmatch_reformat <- function(sub_df, J){
  ntreat <- sum(sub_df$t)
  ncontrol <- sum(1-sub_df$t)
  class <- sub_df$subclass[1]
  
  new_df <- sub_df %>% dplyr::select(y, t, subclass)
  
  # flpi treat and control if ntreat > ncontrol
  if (ntreat > ncontrol){
    flip <- TRUE
    new_df <- arrange(new_df, t)
  }
  else{
    flip <- FALSE
    new_df <- arrange(new_df, -t)
  }
  
  new_df$flip <- flip
  
  #append extra rows
  nextra <- J - ntreat - ncontrol
  
  if (nextra > 0){
    extra_rows <- data.frame(y = rep(NA, nextra), t = NA, subclass = class, flip = flip)
    new_df <- rbind(new_df, extra_rows)
  }
  
  return(new_df)
}