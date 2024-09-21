library(ggplot2)
library(patchwork)
library(tidysdm) # for split violin plots

plot_themes <- function(...) {
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(color = "black", size = 17, margin = margin(0,0,5,0,"mm")),
    legend.text = element_text(color = "black", size = 15),
    ...)
}

violin_plot <- function(plot_data, cov_name) {
  ggplot() +
    geom_split_violin(data = plot_data, aes(x = "", y = !!as.name(cov_name), fill = presabs)) +
    plot_themes(axis.title.x = element_blank())
}

plot_cov_distributions <- function(plot_data, cov_names) {
  cov_plots <- purrr::map(cov_names, \(x) violin_plot(plot_data, x))
  final_plot <- cov_plots[[1]]
  for (i in 2:length(cov_plots)) {
    final_plot <- final_plot + cov_plots[[i]]
  }
  final_plot <- final_plot + plot_layout(guides = "collect")
  return(final_plot)
}

plot_response <- function(cov_name, plot_data, xlabel) {
  ggplot() +
    geom_line(data = plot_data, aes(x = !!as.name(cov_name), y = response)) +
    lims(y = c(0,1))
}

#' Get a tibble for model predictions
#'
#' @param cov_labels List of names for all the covariates in the fitted model
#' @param model Model object
#' @param covname String of the name of the covariate to vary
#' @param fit_data Tibble containing the data used to fit the model
#'
#' @return List containing the gam and lmer objects returned from 
#' gamm4::gamm4
#'
#' @export
get_prediction <- function(cov_labels, model, covname, fit_data) {
  fixed_covs <- cov_labels[cov_labels != covname]
  min_val <- min(fit_data[,covname], na.rm = TRUE)
  max_val <- max(fit_data[,covname], na.rm = TRUE)
  predict_data <- fit_data %>% 
    dplyr::select(all_of(cov_labels)) %>% 
    mutate(across(all_of(fixed_covs), \(x) mean(x, na.rm = T))) %>%
    slice_head(n = 100) %>%
    mutate("{covname}" := seq(min_val, max_val, length.out = 100))
  
  prediction <- predict(model, newdata = predict_data, type = "response") 
  predict_data <- predict_data %>%
    mutate(response = prediction)
  
  return(predict_data)
}

create_pa_response_plots <- function(model, fit_data) {
  cov_names <-  c("depth", "slope", "ice", "cdw", "dshelf")
  covs <- tibble(cov = cov_names) %>%
    mutate(prediction = purrr::map(.$cov, \(x) get_prediction(cov_names, model, x, fit_data)))
  
  plots <- covs %>% mutate(plot = purrr::pmap(., function(cov, prediction, ...) {
    plot_response(cov, prediction)
  })) %>% 
    pull(plot)
  return(plots)
}