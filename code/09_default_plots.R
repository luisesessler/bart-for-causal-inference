results_df <- read.csv(paste0(PATH_RESULTS, "aggregated_metrics_per_dgp.csv"))

results_df$model <- case_when(
  tolower(results_df$model) == "bart_s-learner" ~ "BART/S-learner",
  tolower(results_df$model) == "bart_t-learner" ~ "BART/T-learner",
  tolower(results_df$model) == "bart_ps-glm" ~ "BART/PS-GLM",
  tolower(results_df$model) == "bart_ps-bart" ~ "BART/PS-BART",
  tolower(results_df$model) == "causal_forest" ~ "Causal forest",
  tolower(results_df$model) == "aipw" ~ "AIPW",
  TRUE ~ results_df$model
)

results_df$model <- factor(results_df$model, levels = c(
  "BART/S-learner",
  "BART/T-learner",
  "BART/PS-GLM",
  "BART/PS-BART",
  "Causal forest",
  "AIPW"
))


results_df$dgp <- results_df$dgp %>% str_remove("dgp")


# Boxplots

boxplot_rmse <- ggplot(results_df, aes(x = model, y = rmse_rel)) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.26)) +
  labs(
    x = "Model",
    y = "Relative RMSE"
  ) +
  theme_minimal(base_size = 18)


boxplot_coverage <- ggplot(results_df, aes(x = model, y = coverage)) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0.6, 1)) +
  labs(
    x = "Model",
    y = "Coverage"
  ) +
  theme_minimal(base_size = 18)


boxplot_ci <- ggplot(results_df, aes(x = model, y = ci_rel)) + 
  geom_boxplot() +
  labs(
    x = "Model",
    y = "Relative interval length"
  ) +
  theme_minimal(base_size = 18)


# Scatterplots


scatterplot_rmse <- ggplot(results_df %>% filter(coverage > 0.1), aes(x = dgp, y = rmse_rel,
                                              col = model)) +
  geom_beeswarm(size = 2) +
  labs(x = "DGP", y = "Relative RMSE", col = "Model")+
  theme(
    text = element_text(size = 16),
    legend.position = "top"
  ) +
  scale_color_brewer(palette="Accent")



scatterplot_coverage <- ggplot(results_df %>% filter(coverage > 0.1), aes(x = dgp, y = coverage,
                                              col = model)) +
  geom_beeswarm(size = 2) +
  labs(x = "DGP", y = "Coverage", col = "Model")+
  theme(
    text = element_text(size = 16),
    legend.position = "top"
  ) +
  scale_color_brewer(palette="Accent")



scatterplot_ci <- ggplot(results_df, aes(x = dgp, y = ci_rel,
                   col = model)) +
  geom_beeswarm(size = 2) +
  labs(x = "DGP", y = "Relative interval length", col = "Model")+
  theme(
    text = element_text(size = 16),
    legend.position = "top"
  ) +
  scale_color_brewer(palette="Accent")



plots <- list(
  boxplot_rmse = boxplot_rmse,
  boxplot_coverage = boxplot_coverage,
  boxplot_ci = boxplot_ci,
  scatterplot_rmse = scatterplot_rmse,
  scatterplot_coverage = scatterplot_coverage,
  scatterplot_ci = scatterplot_ci
)

dir.create(PATH_PLOTS, recursive = TRUE, showWarnings = FALSE)

lapply(names(plots), function(name) {
  ggsave(
    filename = paste0(PATH_PLOTS, "/", name, ".png"),
    plot = plots[[name]],
    width = 8, height = 5, dpi = 300, scale = 1.5)
})