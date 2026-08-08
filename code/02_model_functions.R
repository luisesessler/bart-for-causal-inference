# Function for the S-learner models (set method_trt to either "none", "glm" or "bart")
bart_s_learner <- function(y, z, X, true_ate, method_trt = "none", ps_as_covariate = FALSE, 
                           k, n_trees, nu, q){
  
  # fit a single BART model
  bart_fit <- bartc(response = y, treatment = z, confounders = X,
                    method.rsp = "bart", method.trt = method_trt, 
                    estimand = "ate", keepTrees = TRUE, 
                    p.scoreAsCovariate = ps_as_covariate,
                    n.trees = n_trees, k = k,
                    sigdf = nu, sigquant = q)
  

  # Predict outcomes for each observation in treatment and in control group
  predictions_treated <- predict(bart_fit, X, type = "mu.1") # predict outcome under treatment
  predictions_control <- predict(bart_fit, X, type = "mu.0") # predict outcome under control
  predictions_ite <- predictions_treated - predictions_control 
  results_run <- get_metrics_bart(predictions_ite, true_ate)
  return(results_run)
}

# Function for the T-learner
bart_t_learner <- function(y, z, X, true_ate, n_trees, k, nu, q){
  X_mm <- model.matrix(~ . - 1, data = X)
  X_df <- data.frame(X_mm)
  data_mm <- data.frame(X_mm, y)
  
  # Split data into treatment and control group
  treated <- data_mm[z == 1,]
  control <- data_mm[z == 0,]
  
  # define outcome model
  formula_mm <- colnames(X_mm) %>% paste0(collapse = " + ")
  bart_formula <- as.formula(paste0("y ~ ", formula_mm))
  
  # fit two separate BART models
  bart_fit_treated <- bart2(formula = bart_formula, data = treated, n.trees = n_trees, k = k, sigdf = nu, sigquant = q, keepTrees = TRUE, combineChains = TRUE, n.chains = 10)
  bart_fit_control <- bart2(formula = bart_formula, data = control, n.trees = n_trees, k = k, sigdf = nu, sigquant = q, keepTrees = TRUE, combineChains = TRUE, n.chains = 10)
  
  # Prediction for all units in both models
  # "ev" stands for expected value
  predictions_treated <- predict(bart_fit_treated, newdata = X_df, type = "ev") 
  predictions_control <- predict(bart_fit_control, newdata = X_df, type = "ev")
  
  # Calcuate ITEs
  predictions_ite <- predictions_treated - predictions_control
  
  results_run <- get_metrics_bart(predictions_ite, true_ate)
  return(results_run)
}


estimate_causal_forest <- function(y, z, X, true_ate, min_node_size, sample_fraction, mtry){
  X_mm <- model.matrix(~ . - 1, data = X)
  fit_causal_forest <- causal_forest(X_mm, y, z, min.node.size = min_node_size, sample.fraction = sample_fraction, mtry =  mtry)
  
  estimates <- average_treatment_effect(fit_causal_forest)
  ate_estimate <- estimates[1] %>% unname() 
  se <- estimates[2] %>% unname()
  
  results_run <- get_metrics_cf(ate_estimate, se, true_ate)
  return(results_run)
}



aipw <- function(y, z, X, true_ate) {
  
  data <- data.frame(
    y = y,
    z = z,
    X
  )
  
  # Create AIPW object
  fit <- AIPW$new(
    Y = data$y,
    A = data$z,
    W = data[, !(names(data) %in% c("y", "z"))],
    Q.SL.library = "SL.glm", 
    g.SL.library = "SL.glm",
    k_split = 1 
  )
  
  fit$fit()
  
  result <- fit$summary(g.bound = 0.1) # 
  
  ate_est <- fit$result["Mean Difference", "Estimate"]
  lower_ci  <- fit$result["Mean Difference", "95% LCL"]
  upper_ci  <- fit$result["Mean Difference", "95% UCL"]
  
  
  get_metrics_aipw(
    ate_estimate = ate_est,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    ate_true = true_ate
  )
}


# for tuning (BART propensity score model)
aipw_bart <- function(y, z, X, ate_true, g.bound = 0.1) {
  
  X_mm <- model.matrix(~ . - 1, data = X)
  treatment_mm <- data.frame(X_mm, z)
  
  n  <- length(y)
  df <- data.frame(y = y, z = z, X)
  
  # outcome model
  fit <- lm(y ~ z + ., data = df)
  mu1 <- predict(fit, newdata = transform(df, z = 1))
  mu0 <- predict(fit, newdata = transform(df, z = 0))
  
  # treatment model
  g_fit <- bart2(z ~ ., data = treatment_mm, keepTrees = TRUE)
  pred  <- predict(g_fit, newdata = treatment_mm)
  p_score <- colMeans(pred)
  
  # propensity score trucation
  p_score <- pmin(pmax(p_score, g.bound), 1 - g.bound)
  
  ## efficient influence functions
  eif1 <- (z / p_score) * (y - mu1) + mu1
  eif0 <- ((1 - z) / (1 - p_score)) * (y - mu0) + mu0
  
  est <- mean(eif1 - eif0)
  se  <- sd(eif1 - eif0) / sqrt(n)
  ci  <- est + c(-1, 1) * 1.96 * se
  
  get_metrics_aipw(ate_estimate = est, lower_ci = ci[1], upper_ci = ci[2], ate_true = ate_true)
}


# Function for returning the metrics for the BART models
get_metrics_bart <- function(ite_matrix, true_ate){ 
  # Mean ITE for each observation
  ite_estimates <- colMeans(ite_matrix) 
  ate_draws <- rowMeans(ite_matrix)
  
  # ATE Bias
  ate_estimate <- mean(ite_estimates)
  ate_bias <- ate_estimate - true_ate
  
  # Credibility interval
  credible_interval <- quantile(ate_draws, probs = c(0.025, 0.975))
  
  # Interval length
  ci_length <- unname(credible_interval[2] - credible_interval[1])
  
  # Coverage
  is_covered <- ifelse(credible_interval[1] <= true_ate &
                         true_ate <= credible_interval[2], 1, 0) %>% unname()
  
  
  metrics <- c(
    ate_estimate = ate_estimate,
    ate_bias = ate_bias,
    ci_length = ci_length,
    is_covered = is_covered
  )
  return(metrics)
}

# Function for returning the metrics for the causal forest
get_metrics_cf <- function(ate_estimate, se_ate, ate_true){
  # Bias
  ate_bias <- ate_estimate - ate_true
  
  # Confidence intervals
  ci_lower <- ate_estimate - 1.96 * se_ate
  ci_upper <- ate_estimate + 1.96 * se_ate
  
  # CI length
  ci_length <- ci_upper - ci_lower
  
  # Coverage
  is_covered <- ifelse(ci_lower <= ate_true && ate_true <= ci_upper, 1, 0)
  
  # Save results
  metrics <- c(
    ate_estimate = ate_estimate,
    ate_bias = ate_bias,
    ci_length = ci_length,
    is_covered = is_covered
  )
  return(metrics)
}


# Function for returning the metrics for AIPW
get_metrics_aipw <- function(ate_estimate, lower_ci, upper_ci, ate_true){
  # Bias
  ate_bias <- ate_estimate - ate_true
  
  # Confidence intervals
  ci_lower <- lower_ci
  ci_upper <- upper_ci
  
  # CI length
  ci_length <- ci_upper - ci_lower
  
  # Coverage
  is_covered <- ifelse(ci_lower <= ate_true && ate_true <= ci_upper, 1, 0)
  
  metrics <- c(
    ate_estimate = ate_estimate,
    ate_bias = ate_bias,
    ci_length = ci_length,
    is_covered = is_covered
  )
  return(metrics)
}
