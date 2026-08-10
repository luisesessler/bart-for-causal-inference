tuning_dgps <- c(30, 40, 60, 63)
for (dgp_index in tuning_dgps) {
  
  # Read tuning results
  results_df <- read.csv(
    paste0(PATH_RESULTS,
           "tuning/best hyperparameters/best_hyperparameters_vs_default_dgp", 
           dgp_index,
           ".csv"))
  
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
  
  
  
    
    scatterplot_rmse <- ggplot(results_df, aes(x = type, y = rmse_rel,
                                                                          col = model)) +
      geom_beeswarm(size = 2, cex = 2) +
      labs(x = "Setting", y = "Relative RMSE", col = "Model")+
      theme(
        text = element_text(size = 16),
        legend.position = "top"
      ) +
      scale_color_brewer(palette="Accent")
    
    name <- paste0("scatterplot_tuning_", dgp_index)
    
    ggsave(
      filename = paste0(PATH_PLOTS, "/", name, ".png"),
      scatterplot_rmse,
      width = 4, height = 5, dpi = 300, scale = 1.5)
}