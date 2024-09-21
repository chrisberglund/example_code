library(tidyverse)
library(gamm4)
library(dismo)
source("4_figures.R")

#' Wrapper function for the generalized additive mixed model 
#'
#' @param fit_data Tibble containing the data to fit the model to
#' @param weight The weights for each row
#'
#' @return List containing the gam and lmer objects returned from 
#' gamm4::gamm4
#'
#' @export
pa_gamm <- function(fit_data, weight, ...) {
  gamm4(presabs ~ s(depth, k = 5, bs = "ts") + 
               s(cdw, k = 5, bs = "ts") + 
               s(ice, k = 5, bs = "ts") + 
               s(dshelf, k = 5, bs = "ts") + 
               s(slope, k = 5, bs = "ts"),
             random = ~ (1 | id),
             weights = weight,
             family = binomial,
             data = fit_data,
       ...)
}

#' Wrapper function for the boosted regression tree 
#'
#' @param fit_data Tibble containing the data to fit the model to
#' @param weight The weights for each row
#'
#' @return 
#'
#' @export
pa_brt <- function(fit_data, ...) {
  input <- fit_data %>% dplyr::select(presabs, depth, slope, ice, cdw, dshelf)
  # use hyperparameters determined from cross validation
  gbm.step(data = as.data.frame(input), gbm.x = c(2,3,4,5,6), gbm.y = 1, 
           tree.complexity = 2, 
           learning.rate = 0.05,
           bag.fraction = 0.75, 
           family = "bernoulli",
           ...)
}

gamm_summary <- function(model) {
  fit_data <- model$gam$model
  null_model <- gamm4(presabs ~ 1,
                      random = ~ (1 |id),
                      weights = unlist(fit_data["(weights)"]),
                      family = binomial,
                      data = fit_data)
  null_deviance <- deviance(null_modl$mer)
  resid_deviance <- deviance(model$mer)
  deviance_explained <- (null_deviance - resid_deviance) / null_deviance
  return(deviance_explained)
}

#' Get bootstraped confidence intervals for boosted regression tree predictions
#'
#' @param model gbm Object returned from the fitted model
#' @param predict_df tibble Contains values to be used for predictions
#' @param input_data tibble Contains values used for fitting the original model (columns must be in same order)
#'
#' @return tibble Contains predictions and associated confidence intervals on the natural scale
#'
#' @export
bootstrap_brt_ci <- function(model, predict_df, input_data, nreplicates = 100) {
  boot_predicts <- matrix(0, ncol = nreplicates, nrow = nrow(predict_df))
  for (i in 1:nreplicates) {
    rep_data <- input_data %>% slice_sample(n =nrow(input_data), replace = TRUE)
    
    weight <- rep_data$weight
    fit <- gbm.fixed(data = rep_data,
                     gbm.x = c(2,3,4,5,6), 
                     gbm.y = 1,
                     n.trees = model$n.trees,
                     tree.complexity = model$interaction.depth,
                     learning.rate = model$shrinkage,
                     bag.fraction = model$bag.fraction)
    boot_predicts[,i] <- predict(fit, newdata = predict_df, type = "response")
  }
  ci = apply(boot_predicts, 1, quantile,c(0.025, 0.975))
  predict_df <- predict_df %>%
    dplyr::mutate(response = model$family$linkinv(predict(model, newdata=., type = "response")),
                  lwr = model$family$linkinv(ci[1,]),
                  upr = model$family$linkinv(ci[2,]))
}

#' Get bootstraped confidence intervals for boosted regression tree predictions
#'
#' @param model gbm Object returned from the fitted model
#' @param predict_df tibble Contains values to be used for predictions
#' @param input_data tibble Contains values used for fitting the original model (columns must be in same order)
#'
#' @return tibble Contains predictions and associated confidence intervals on the natural scale
#'
#' @export
bootstrap_brt_ci <- function(model, predict_df, input_data, nreplicates = 100) {
  boot_predicts <- matrix(0, ncol = nreplicates, nrow = nrow(predict_df))
  for (i in 1:nreplicates) {
    rep_data <- input_data %>% slice_sample(n =nrow(input_data), replace = TRUE)
    
    weight <- rep_data$weight
    fit <- gbm.fixed(data = rep_data,
                     gbm.x = c(2,3,4,5,6), 
                     gbm.y = 1,
                     n.trees = model$n.trees,
                     tree.complexity = model$interaction.depth,
                     learning.rate = model$shrinkage,
                     bag.fraction = model$bag.fraction)
    boot_predicts[,i] <- predict(fit, newdata = predict_df, type = "response")
  }
  ci = apply(boot_predicts, 1, quantile,c(0.025, 0.975))
  predict_df <- predict_df %>%
    dplyr::mutate(response = predict(model, newdata=., type = "response"),
                  lwr = ci[1,],
                  upr = ci[2,])
}


#' Separates data into a training set and a test set
#'
#' @param pa_data tibble The dataset to split
#' @param ratio numeric The percentage of rows to be assigned to the training set
#'
#' @return tibble pa_data with an added boolean column where TRUE is for the training set
#'
#' @export
separate_sets <- function(pa_data, ratio) {
  pa_data <- pa_data %>% 
    group_by(id, datetime) %>% 
    mutate(step = cur_group_id())
  steps <- sample(1:max(pa_data$step), ratio * max(pa_data$step), replace = FALSE)
  pa_data <- mutate(pa_data, fit = step %in% steps)
  return(pa_data)
}

### load data ###

eseal_data <- readRDS("data/eseal_pa_data.RDS")
wseal_data <- readRDS("data/wseal_pa_data.RDS") 

eseal_gamm <- pa_gamm(eseal_data, eseal_data$weight)
wseal_gamm <- pa_gamm(wseal_data, wseal_data$weight)

eseal_brt <- pa_brt(eseal_data)
wseal_brt <- pa_brt(wseal_brt)
