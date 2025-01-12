# Load packages ----
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(
  cmdstanr,
  posterior,
  loo,
  qs
)

source("Scripts/load_data.R")

loo_list <- list()
Model_M <- cmdstan_model(file.path("Scripts/models", "mLP_1.stan"), force_recompile = T) # For moment matching
Model_PH <- cmdstan_model(file.path("Scripts/models", "mPH_1.stan"), force_recompile = T)
Model_H <- cmdstan_model(file.path("Scripts/models", "model_hybrid.stan"), force_recompile = T)


# PH model ----
# Dataset VMAC: 
if (!file.exists("Output/results/model_PH_VMAC.qs")) {
  fit <- Model_PH$sample(
    data = dStan[["VMAC"]],
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1500,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    init = generate_initial_values(J = dStan$VMAC$J, model = "2")
  )
  fit$draws()
  qsave(fit, "Output/results/model_PH_VMAC.qs")
  
  if (!file.exists("Output/results/loo_estimates.rds")) {
    loo_list[["PH_VMAC"]] <- fit$loo(moment_match = T) # some bad k_hat values (close to 0.7 threshold)
  }
}

# Dataset UMAC
if (!file.exists("Output/results/model_PH_UMAC.qs")) {
  fit <- Model_PH$sample(
    data = dStan[["UMAC"]],
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1500,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    init = generate_initial_values(J = dStan$UMAC$J, model = "2")
  )
  fit$draws()
  qsave(fit, "Output/results/model_PH_UMAC.qs")
  
  if (!file.exists("Output/results/loo_estimates.rds")) {
    loo_list[["PH_UMAC"]] <- fit$loo()
  }
}

# Mackintosh model ----
# Dataset VMAC: 
if (!file.exists("Output/results/model_M_VMAC.qs")) {
  fit <- Model_M$sample(
    data = dStan[["VMAC"]],
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1500,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    init = generate_initial_values(J = dStan$VMAC$J)
  )
  fit$draws()
  qsave(fit, "Output/results/model_M_VMAC.qs")
  
  if (!file.exists("Output/results/loo_estimates.rds")) {
    loo_list[["M_VMAC"]] <- fit$loo(moment_match = T) # some bad k_hat values (close to 0.7 threshold)
  }
}

# Dataset UMAC
if (!file.exists("Output/results/model_M_UMAC.qs")) {
  fit <- Model_M$sample(
    data = dStan[["UMAC"]],
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1500,
    iter_sampling = 2000,
    adapt_delta = 0.95,
    init = generate_initial_values(J = dStan$UMAC$J)
  )
  fit$draws()
  qsave(fit, "Output/results/model_M_UMAC.qs")
  
  if (!file.exists("Output/results/loo_estimates.rds")) {
    loo_list[["M_UMAC"]] <- fit$loo()
  }
}

# Save PSIS-LOO: ----
if(!file.exists("Output/results/loo_estimates.rds")) saveRDS(loo_list, "Output/results/loo_estimates.rds")

# PPC: ----
if(!file.exists("Output/figures/dPlots.rds")) {
  predList <- list()
  
  # PH UMAC
  fit <- qread("Output/results/model_PH_UMAC.qs")
  df <- get_PPC(fit)
  
  predList[["PH_UMAC"]] <- df %>%
    mutate(Trial = rep(1:512, length.out = nrow(.)),
           Singleton = rep(c("Low", "High"), each = nrow(.)/2),
           Epoch = create_epochs(Trial, 32)
    ) %>%
    dplyr::summarise(P = median(Mean),
                     Low = mean(Low),
                     High = mean(High),
                     .by = c(Epoch, Singleton)) %>%
    mutate(Model = "Pearce-Hall",
           Data = "Chow et al. (2024)")
  
  rm(fit)
  
  # M UMAC
  fit <- qread("Output/results/model_M_UMAC.qs")
  df <- get_PPC(fit)
  
  predList[["M_UMAC"]] <- df %>%
    mutate(Trial = rep(1:512, length.out = nrow(.)),
           Singleton = rep(c("Low", "High"), each = nrow(.)/2),
           Epoch = create_epochs(Trial, 32)
    ) %>%
    dplyr::summarise(P = median(Mean),
                     Low = mean(Low),
                     High = mean(High),
                     .by = c(Epoch, Singleton)) %>%
    mutate(Model = "Mackintosh",
           Data = "Chow et al. (2024)")
  
  # PH VMAC
  fit <- qread("Output/results/model_PH_VMAC.qs")
  df <- get_PPC(fit, group_size = 384)
  
  predList[["PH_VMAC"]] <- df %>%
    mutate(Trial = rep(1:384, length.out = nrow(.)),
           Singleton = rep(c("High", "Low"), each = nrow(.)/2),
           Epoch = create_epochs(Trial, 32)
    ) %>%
    dplyr::summarise(P = median(Mean),
                     Low = mean(Low),
                     High = mean(High),
                     .by = c(Epoch, Singleton)) %>%
    mutate(Model = "Pearce-Hall",
           Data = "Le et al. (2024)")
  
  # M VMAC
  fit <- qread("Output/results/model_M_VMAC.qs")
  df <- get_PPC(fit, group_size = 384)
  
  predList[["M_VMAC"]] <- df %>%
    mutate(Trial = rep(1:384, length.out = nrow(.)),
           Singleton = rep(c("High", "Low"), each = nrow(.)/2),
           Epoch = create_epochs(Trial, 32)
    ) %>%
    dplyr::summarise(P = median(Mean),
                     Low = mean(Low),
                     High = mean(High),
                     .by = c(Epoch, Singleton)) %>%
    mutate(Model = "Mackintosh",
           Data = "Le et al. (2024)")
  
  rm(fit, df) # rm unused objects
  
  # Save data
  saveRDS(list(do.call(rbind, predList),
               do.call(rbind, dPlots)), "Output/figures/dPlots.rds")
}