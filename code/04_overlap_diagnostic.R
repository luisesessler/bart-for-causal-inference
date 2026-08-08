dgps <- read.csv(PATH_DGP_OVERVIEW)
all_datasets <- read.csv(PATH_FILE_OVERVIEW)

set.seed(213)


overlap_bart <- data.frame(dgp = character(),
           iteration = numeric(),
           n_discarded = numeric())

for (i in 1:nrow(dgps)){
  dgp_index <- dgps$dgp[i]
  datasets_dgp <- all_datasets %>% filter(DGPid == dgp_index)

  for (j in 1:5){
    file_name <- datasets_dgp[j, "filename"]
    
  data <- read.csv(paste0(PATH_DATASETS, file_name, ".csv"))
  y <- data$Y
  z <- data$A
  X <- data %>% select(starts_with("V"))
  
  # Fit BART model
  bart_fit <- bartc(
    response = y,
    treatment = z,
    confounders = X,
    method.rsp = "bart",
    method.trt = "none",
    estimand = "ate",
    keepTrees = TRUE,
    p.scoreAsCovariate = FALSE
  )
  
  
  # Posterior draws
  mu1_draws <- predict(bart_fit, X, type = "mu.1")
  mu0_draws <- predict(bart_fit, X, type = "mu.0")
  
  # Posterior SDs
  s1 <- apply(mu1_draws, 2, sd)  # SD(mu.1)
  s0 <- apply(mu0_draws, 2, sd)  # SD(mu.0)
  
  # --- Max rule ---
  # Factual SDs within groups
  s1_treated <- s1[z == 1]  # observed for treated
  s0_control <- s0[z == 0]  # observed for control
  
  # Thresholds: max + buffer (SD of SDs)
  threshold_treated <- max(s1_treated) + sd(s1_treated)
  threshold_control <- max(s0_control) + sd(s0_control)
  
  # Initialize
  discard_treated <- rep(FALSE, length(z))
  discard_control <- rep(FALSE, length(z))
  
  # For treated units: use counterfactual SD = s0
  discard_treated[z == 1] <- s0[z == 1] > threshold_treated
  
  # For control units: use counterfactual SD = s1
  discard_control[z == 0] <- s1[z == 0] > threshold_control
  
  # Combine (units to discard)
  discard <- discard_treated | discard_control
  
  overlap_bart <- rbind(overlap_bart, data.frame(
    dgp         = dgp_index,
    iteration   = j,
    n_discarded = sum(discard)
  ))
  }
}

write.csv(overlap_bart, paste0(PATH_RESULTS, "diagnostics/overlap_diagnostic.csv"), row.names = FALSE)


# Summary
overlap_summary <- overlap_bart %>%
  group_by(dgp) %>%
  summarise(mean_n_discarded = mean(n_discarded)) %>%
  left_join(dgps %>% select(dgp, n_observation), by = "dgp") %>%
  mutate(discarded_ratio = mean_n_discarded / n_observation)

write.csv(overlap_summary, paste0(PATH_RESULTS, "diagnostics/overlap_diagnostic_summary.csv"), row.names = FALSE)

