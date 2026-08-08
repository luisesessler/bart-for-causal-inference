dgps <- read.csv(PATH_DGP_OVERVIEW)

true_ate_lookup <- dgps %>%
  transmute(
    dgp = paste0("dgp", dgp),
    trueATE
  )

results_df <- read.csv(paste0(PATH_RESULTS, "results_default_all.csv"))


# Summaries per DGP and model

result_aggregated <- results_df %>%
  group_by(dgp, model, metric) %>%
  summarise(
    value = case_when(
      first(metric) == "ate_bias" ~ sqrt(mean(value^2)),
      TRUE                        ~ mean(value)
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
    rmse_rel = ate_bias / abs(trueATE),
    coverage,
    ci_rel = ci_length / abs(trueATE)
  )


write.csv(result_aggregated, paste0(PATH_RESULTS, "aggregated_metrics_per_dgp.csv"), row.names = FALSE)


# Summaries per model

result_per_model <- result_aggregated %>%
  group_by(model) %>%
  summarise(
    rmse_rel = mean(rmse_rel),
    coverage = mean(coverage),
    ci_rel = mean(ci_rel),
    .groups = "drop"
  )

write.csv(result_per_model, paste0(PATH_RESULTS, "aggregated_metrics_per_model.csv"), row.names = FALSE)
