dgps <- read.csv(PATH_DGP_OVERVIEW)
all_datasets <- read.csv(PATH_FILE_OVERVIEW)

set.seed(213)

# DGPs and dataset subset used for tuning
tuning_dgps <- c(30, 40, 60, 63)
n_datasets_tuning  <- 20  # only the first 20 datasets per DGP are used for tuning

# BART hyperparameter grid

sigma_params <- list(c(3, 0.90), c(3, 0.99), c(10, 0.75))
k_vals  <- c(1, 2, 3, 5)
n_trees_vals <- c(75, 100, 200)

hyper_bart <- expand.grid(sigma_idx = 1:3, k = k_vals, n_trees = n_trees_vals)
hyper_bart$nu <- sapply(hyper_bart$sigma_idx, function(idx) sigma_params[[idx]][1])
hyper_bart$q  <- sapply(hyper_bart$sigma_idx, function(idx) sigma_params[[idx]][2])
hyper_bart$sigma_idx <- NULL


# ---- BART meta-learner variants that share the hyper_bart grid -------------
# (T-learner is handled separately below, since bart_t_learner() has a
# different function signature and does not take method_trt/ps_as_covariate)
bart_s_learners <- list(
  list(model_name = "bart_s-learner", method_trt = "none", ps_as_covariate = FALSE),
  list(model_name = "bart_ps-glm", method_trt = "glm",  ps_as_covariate = TRUE),
  list(model_name = "bart_ps-bart", method_trt = "bart", ps_as_covariate = TRUE)
)

# Causal forest hyperparameter grid 

min_node_size_vals  <- c(5, 10, 20, 50)
sample_fraction_vals <- c(0.3, 0.5)

# Iterate through all tuning DGPs
for (dgp_index in tuning_dgps) {
  
  # Initialize results data frame
  df_tuning_results <- data.frame(
    dgp = character(),
    metric= character(),
    model = character(),
    iteration = numeric(),
    k = numeric(),
    n_trees = numeric(),
    nu = numeric(),
    q = numeric(),
    min_node_size = numeric(),
    sample_fraction = numeric(),
    mtry  = numeric(),
    weights_method = character(),
    value = numeric()
  )
  
  
  
  dgp_row  <- which(dgps$dgp == dgp_index)
  true_ate <- dgps$trueATE[dgp_row]
  
  datasets_dgp <- all_datasets %>% filter(DGPid == dgp_index)
  
  # Iterate through 20 datasets
  for (j in 1:n_datasets_tuning) {
    
    file_name <- datasets_dgp[j, "filename"]
    data <- read.csv(paste0(PATH_DATASETS, file_name, ".csv"))
    y <- data$Y
    z <- data$A
    X <- data %>% select(starts_with("V"))
    
    # Dynamically set mtry
    n_cov <- ncol(X)
    hyper_cf <- expand.grid(
      min_node_size = min_node_size_vals,
      sample_fraction = sample_fraction_vals,
      mtry = c(
        ceiling(n_cov / 2),
        ceiling(n_cov / 3),
        min(ceiling(sqrt(n_cov)) + 20, n_cov)
      )
    )
    
    # Causal forest
    for (l in 1:nrow(hyper_cf)) {
      hp <- hyper_cf[l, ]
      
      metrics_run <- estimate_causal_forest(
        y, z, X, true_ate,
        min_node_size = hp$min_node_size,
        sample_fraction = hp$sample_fraction,
        mtry = hp$mtry
      )
      
      if (length(metrics_run) != 4) {
        cat("Malformed metrics_run at dgp:", dgp_index, "j:", j, "l:", l, "\n")
        cat("min_node_size:", hp$min_node_size, "sample_fraction:", hp$sample_fraction, "mtry:", hp$mtry, "\n")
        print(metrics_run)
        stop("metrics_run length != 4 — see printed details above")
      }
      
      df_tuning_results <- rbind(
        df_tuning_results,
        data.frame(
          dgp = rep(paste0("dgp", dgp_index), 4),
          metric = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
          model = rep("causal_forest", 4),
          iteration = rep(j, 4),
          k = NA,
          n_trees = NA,
          nu = NA,
          q = NA,
          min_node_size = rep(hp$min_node_size, 4),
          sample_fraction = rep(hp$sample_fraction, 4),
          mtry = rep(hp$mtry, 4),
          weights_method = rep(NA, 4),
          value  = metrics_run
        )
      )
    }

    # S-learners
    for (variant in bart_s_learners) {
      for (l in 1:nrow(hyper_bart)) {
        hp <- hyper_bart[l, ]
        
        metrics_run <- bart_s_learner(
          y, z, X, true_ate = true_ate,
          method_trt      = variant$method_trt,
          ps_as_covariate = variant$ps_as_covariate,
          k = hp$k, n_trees = hp$n_trees, nu = hp$nu, q = hp$q
        )
        
        df_tuning_results <- rbind(
          df_tuning_results,
          data.frame(
            dgp = rep(paste0("dgp", dgp_index), 4),
            metric = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
            model = rep(variant$model_name, 4),
            iteration = rep(j, 4),
            k  = rep(hp$k, 4),
            n_trees = rep(hp$n_trees, 4),
            nu = rep(hp$nu, 4),
            q = rep(hp$q, 4),
            min_node_size = NA,
            sample_fraction = NA,
            mtry  = NA,
            weights_method = rep(NA, 4),
            value = metrics_run
          )
        )
      }
    }
    
    # BART T-learner
    for (l in 1:nrow(hyper_bart)) {
      hp <- hyper_bart[l, ]
      
      metrics_run <- bart_t_learner(
        y, z, X, true_ate = true_ate,
        k = hp$k, n_trees = hp$n_trees, nu = hp$nu, q = hp$q
      )
      
      df_tuning_results <- rbind(
        df_tuning_results,
        data.frame(
          dgp = rep(paste0("dgp", dgp_index), 4),
          metric  = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
          model = rep("bart_t-learner", 4),
          iteration = rep(j, 4),
          k = rep(hp$k, 4),
          n_trees = rep(hp$n_trees, 4),
          nu = rep(hp$nu, 4),
          q = rep(hp$q, 4),
          min_node_size = NA,
          sample_fraction = NA,
          mtry = NA,
          weights_method  = rep(NA, 4),
          value = metrics_run
        )
      )
    }
    
    # AIPW
    metrics_iteration_aipw <- aipw_bart(y = y, z = z, X = X, ate_true = true_ate)
    
    df_tuning_results <- rbind(
      df_tuning_results,
      data.frame(
        dgp = rep(paste0("dgp", dgp_index), 4),
        metric = c("ate_estimate", "ate_bias", "ci_length", "coverage"),
        model  = rep("AIPW", 4),
        iteration = rep(j, 4),
        k = NA,
        n_trees = NA,
        nu = NA,
        q = NA,
        min_node_size = NA,
        sample_fraction = NA,
        mtry = NA,
        weights_method = rep("bart", 4),
        value = metrics_iteration_aipw
      )
    )
    print(paste("dgp:", dgp_index, "- finished dataset", j, "of", n_datasets_tuning))
  }
  
  # Save after each DGP finishes tuning, so results for completed DGPs are
  # available even if a later DGP is not yet run
  write.csv(df_tuning_results,
            file = paste0(PATH_RESULTS, "tuning/tuning_results_all_models_dgp", dgp_index, ".csv"),
            row.names = FALSE)
  
}