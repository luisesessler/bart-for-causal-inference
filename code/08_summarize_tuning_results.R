dgps <- read.csv(PATH_DGP_OVERVIEW)

true_ate_lookup <- dgps %>%
  transmute(
    dgp = paste0("dgp", dgp),
    trueATE
  )

tuning_dgps <- c(30, 40, 60, 63)
for (dgp_index in tuning_dgps) {
  
  # Read tuning results
  results_df <- read.csv(
    paste0(PATH_RESULTS,
           "tuning/tuning_results_all_models_dgp",
           dgp_index,
           ".csv")
  )
  
  # Aggregate over the 20 datasets
  result_aggregated <- results_df %>%
    group_by(
      dgp,
      model,
      k,
      n_trees,
      nu,
      q,
      min_node_size,
      sample_fraction,
      mtry,
      weights_method,
      metric
    ) %>%
    summarise(
      value = case_when(
        first(metric) == "ate_bias" ~ sqrt(mean(value^2)),
        TRUE ~ mean(value)
      ),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = metric,
      values_from = value
    ) %>%
    left_join(true_ate_lookup, by = "dgp") %>%
    transmute(
      dgp,
      model,
      k,
      n_trees,
      nu,
      q,
      min_node_size,
      sample_fraction,
      mtry,
      weights_method,
      rmse_rel = ate_bias / abs(trueATE),
      coverage,
      ci_rel = ci_length / abs(trueATE)
    )
  
  # Remove irrelevant tuning columns for each model
  result_aggregated <- result_aggregated %>%
    mutate(
      k = ifelse(grepl("^bart", model), k, NA),
      n_trees = ifelse(grepl("^bart", model), n_trees, NA),
      nu = ifelse(grepl("^bart", model), nu, NA),
      q = ifelse(grepl("^bart", model), q, NA),
      min_node_size = ifelse(model == "causal_forest", min_node_size, NA),
      sample_fraction = ifelse(model == "causal_forest", sample_fraction, NA),
      mtry = ifelse(model == "causal_forest", mtry, NA),
      weights_method = ifelse(model == "AIPW", weights_method, NA)
    )
  
  # Order rows
  result_aggregated <- result_aggregated %>%
    arrange(
      model,
      k,
      n_trees,
      nu,
      q,
      min_node_size,
      sample_fraction,
      mtry,
      weights_method
    )
  
  # Save all aggregated results
  write.csv(
    result_aggregated,
    file = paste0(
      PATH_RESULTS,
      "tuning/summaries/tuning_summary_dgp",
      dgp_index,
      ".csv"
    ),
    row.names = FALSE
  )
  
  # Select the best hyperparameters
  best_hyperparameters <- result_aggregated %>%
    group_by(model) %>%
    slice_min(
      order_by = rmse_rel,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  
  # Add default run results for comparison
  # assumes results_default_all.csv has the same structure as the tuning
  # results (dgp, metric, model, iteration, value), and that iteration
  # order matches the same datasets used for tuning
  default_df <- read.csv(paste0(PATH_RESULTS, "results_default_all.csv")) %>%
    filter(dgp == paste0("dgp", dgp_index), iteration <= 20)
  
  default_aggregated <- default_df %>%
    group_by(dgp, model, metric) %>%
    summarise(
      value = case_when(
        first(metric) == "ate_bias" ~ sqrt(mean(value^2)),
        TRUE ~ mean(value)
      ),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = metric, values_from = value) %>%
    left_join(true_ate_lookup, by = "dgp") %>%
    transmute(
      dgp,
      model,
      k = ifelse(grepl("^bart", model), 2, NA),
      n_trees = ifelse(grepl("^bart", model), 75, NA),
      nu = ifelse(grepl("^bart", model), 3, NA),
      q = ifelse(grepl("^bart", model), 0.9, NA),
      min_node_size = ifelse(model == "causal_forest", 5, NA),
      sample_fraction = ifelse(model == "causal_forest", 0.5, NA),
      mtry = ifelse(model == "causal_forest",
                    pmin(ceiling(sqrt(dgps$n_covariates[dgps$dgp == dgp_index])) + 20,
                         dgps$n_covariates[dgps$dgp == dgp_index]),
                    NA),
      weights_method = ifelse(model == "AIPW", "glm", NA), # check this matches your naming for the default propensity model
      rmse_rel = ate_bias / abs(trueATE),
      coverage,
      ci_rel = ci_length / abs(trueATE),
      type = "default"
    )
  
  best_hyperparameters <- best_hyperparameters %>%
    mutate(type = "tuned")
  
  combined_results <- bind_rows(best_hyperparameters, default_aggregated) %>%
    arrange(model, type)
  
  write.csv(
    combined_results,
    file = paste0(
      PATH_RESULTS,
      "tuning/best hyperparameters/best_hyperparameters_vs_default_dgp",
      dgp_index,
      ".csv"
    ),
    row.names = FALSE
  )
}