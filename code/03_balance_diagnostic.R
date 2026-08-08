dgps <- read.csv(PATH_DGP_OVERVIEW)
all_datasets <- read.csv(PATH_FILE_OVERVIEW)

balance_results <- data.frame(dgp = character(),
                              iteration = numeric(),
                              balance = numeric())

for (i in 1:nrow(dgps)) {
  dgp_index <- dgps$dgp[i]
  datasets_dgp <- all_datasets %>% filter(DGPid == dgp_index)
  
  for (j in 1:5) {
    file_name <- datasets_dgp[j, "filename"]
    data <- read.csv(paste0(PATH_DATASETS, file_name, ".csv"))
    z <- data$A
    
    
    balance_results <- rbind(balance_results, data.frame(
      dgp = dgp_index,
      iteration = j,
      balance = mean(z)
      )
    )
      
  }
}

write.csv(balance_results, paste0(PATH_RESULTS, "diagnostics/balance_diagnostic.csv"), row.names = FALSE)

balance_summary <- balance_results %>%
  group_by(dgp) %>%
  summarise(mean_balance = mean(balance)) %>%
  left_join(dgps %>% select(dgp, n_observation), by = "dgp") %>%
  mutate(n_smaller_group = round(pmin(mean_balance, 1 - mean_balance) * n_observation))

write.csv(balance_summary, paste0(PATH_RESULTS, "diagnostics/balance_summary.csv"), row.names = FALSE)
