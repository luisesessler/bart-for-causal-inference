dgps <- read.csv(PATH_DGP_OVERVIEW)
all_datasets <- read.csv(PATH_FILE_OVERVIEW)

set.seed(213)

# Default hyperparameters
# BART defaults 
default_k       <- 2
default_n_trees <- 75
default_nu      <- 3
default_q       <- 0.90

# Causal forest default
# mtry is dynamically set based on the number of covariates, so it will be set within the loop
default_min_node_size  <- 5
default_sample_fraction <- 0.5


# Initialize results dataframe
df_results_default <- data.frame(
  dgp       = character(),
  metric    = character(),
  model     = character(),
  iteration = numeric(),
  value     = numeric()
)

# Iterate through all 16 DGPs
for(i in 1:nrow(dgps)){
  dgp_index <- dgps$dgp[i]
  true_ate <- dgps$trueATE[i]
  
  # Get all file names of each of the datasets of the DGP
  datasets_dgp <- all_datasets %>% filter(DGPid == dgp_index)
  
  # Iterate through all 100 dataset of each DGP
  for (j in 1:nrow(datasets_dgp)){
    file_name <- datasets_dgp[j, "filename"]
    data <- read.csv(paste0(PATH_DATASETS, file_name, ".csv"))
    y <- data$Y
    z <- data$A
    X <- data %>% select(starts_with("V"))
    
    # S-learner
    metrics_iteration_s_learner <- bart_s_learner(y, z, X, true_ate = true_ate, method_trt = "none",
                     k = default_k, n_trees = default_n_trees,
                     nu = default_nu, q = default_q)
    
    df_results_default <- rbind(df_results_default, data.frame(
      dgp       = rep(paste0("dgp", dgp_index), 4),
      metric    = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
      model     = rep("bart_s-learner", 4),
      iteration = rep(j, 4),
      value     = metrics_iteration_s_learner
    ))
    
    
    # T-learner
    metrics_iteration_t_learner <- bart_t_learner(y, z, X, true_ate = true_ate,
                     k = default_k, n_trees = default_n_trees,
                     nu = default_nu, q = default_q)
    
    df_results_default <- rbind(df_results_default, data.frame(
      dgp       = rep(paste0("dgp", dgp_index), 4),
      metric    = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
      model     = rep("bart_t-learner", 4),
      iteration = rep(j, 4),
      value     = metrics_iteration_t_learner
    ))
    
    
    
    # PS-GLM
    metrics_iteration_ps_glm <- bart_s_learner(y, z, X, true_ate = true_ate, method_trt = "glm", ps_as_covariate = TRUE,
                     k = default_k, n_trees = default_n_trees,
                     nu = default_nu, q = default_q)

    
    df_results_default <- rbind(df_results_default, data.frame(
      dgp       = rep(paste0("dgp", dgp_index), 4),
      metric    = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
      model     = rep("bart_ps-glm", 4),
      iteration = rep(j, 4),
      value     = metrics_iteration_ps_glm
    ))
    
    
    # PS-BART
    metrics_iteration_ps_bart <- bart_s_learner(y, z, X, true_ate = true_ate, method_trt = "bart", ps_as_covariate = TRUE,
                     k = default_k, n_trees = default_n_trees,
                     nu = default_nu, q = default_q)
    
    df_results_default <- rbind(df_results_default, data.frame(
      dgp       = rep(paste0("dgp", dgp_index), 4),
      metric    = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
      model     = rep("bart_ps-bart", 4),
      iteration = rep(j, 4),
      value     = metrics_iteration_ps_bart
    ))
    
    
    # Causal forest
    # dynamically set mtry
    p_cov <- ncol(X)
    default_mtry <- min(ceiling(sqrt(p_cov) + 20), p_cov)
    
    metrics_iteration_cf <- estimate_causal_forest(
        y, z, X, true_ate,
        min_node_size   = default_min_node_size,
        sample_fraction = default_sample_fraction,
        mtry            = default_mtry
      )
    
    df_results_default <- rbind(df_results_default, data.frame(
      dgp       = rep(paste0("dgp", dgp_index), 4),
      metric    = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
      model     = rep("causal_forest", 4),
      iteration = rep(j, 4),
      value     = metrics_iteration_cf
    ))
    
    
    # AIPW
    metrics_iteration_aipw <- aipw(y = y, z = z, X = X, true_ate = true_ate)
    
    df_results_default <- rbind(df_results_default, data.frame(
      dgp       = rep(paste0("dgp", dgp_index), 4),
      metric    = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
      model     = rep("AIPW", 4),
      iteration = rep(j, 4),
      value     = metrics_iteration_aipw
    ))
    
  }
}

write.csv(df_results_default,
          file = paste0(PATH_RESULTS, "results_default_all.csv"),
          row.names = FALSE)
