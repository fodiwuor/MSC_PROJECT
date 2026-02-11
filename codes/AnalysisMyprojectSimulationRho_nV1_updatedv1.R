
## =========================================================
##  BLOCK 1: Build performance summary from AllRhoCombined
##  (PRETREND + CORRELATIONS REMOVED as requested)
## =========================================================
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

library(MASS)
library(dplyr)
library(sandwich)
library(lmtest)
library(purrr)
library(ggplot2)

## --- 1. Load your big simulated data set ---
AllRhoCombined <- readRDS("data/AllRhoCombined_diffAR_ZvsY.rds")
str(AllRhoCombined)

## --- 2. True policy effects for each scenario (on log scale) ---
true_B_fun <- function(scen) {
  dplyr::case_when(
    scen == "Small"     ~ -0.0408,
    scen == "Moderate"  ~ -0.3567,
    scen == "Large"     ~ -0.5108,
    TRUE                ~ NA_real_
  )
}

## --- 3. Model settings ---
delta  <- 1e-4                          # avoid log(0)
fam_nb <- MASS::negative.binomial(350)  # same dispersion as simulation

## --- 4. Fit one simulation for a given method (NO pretrend/cor stuff) ---
fit_one_sim <- function(dat, method = c("Trd", "CITS", "CITS_0.4_constant", "CITS_spill_3pct")) {
  method <- match.arg(method)
  L <- floor(nrow(dat)^(1/4))   # HAC lag rule
  
  fit <- tryCatch(
    {
      if (method == "CITS") {
        glm(Y ~ P + offset(log(Z + delta)),
            family = fam_nb, data = dat)
        
      } else if (method == "CITS_0.4_constant") {
        glm(Y ~ P + offset(log(Z_zeropoint4Cor + delta)),
            family = fam_nb, data = dat)
        
      } else if (method == "CITS_spill_3pct") {
        glm(Y ~ P + offset(log(ZsameDirSpil_3pct + delta)),
            family = fam_nb, data = dat)
        
      } else {
        glm(Y ~ t + P + X1,
            family = fam_nb, data = dat)
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, converged = FALSE))
  }
  
  if (!isTRUE(fit$converged)) {
    return(list(estimate = NA_real_, se = NA_real_, p = NA_real_, converged = FALSE))
  }
  
  out <- tryCatch(
    {
      V <- sandwich::NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE)
      
      est <- coef(fit)["P"]
      se  <- sqrt(diag(V)["P"])
      z   <- est / se
      p   <- 2 * pnorm(-abs(z))
      
      list(estimate = est, se = se, p = p, converged = TRUE)
    },
    error = function(e) {
      est <- coef(fit)["P"]
      list(estimate = est, se = NA_real_, p = NA_real_, converged = FALSE)
    }
  )
  
  out
}

## --- 5. Simulation-level estimates (rho, n, scenario, j, method) ---
estimates <- AllRhoCombined %>%
  arrange(rho, n, EstimandScenario, j, t) %>%
  group_by(rho, n, EstimandScenario, j) %>%
  group_modify(~{
    dat <- .x
    
    res_trd     <- fit_one_sim(dat, "Trd")
    res_cits    <- fit_one_sim(dat, "CITS")
    res_cits04  <- fit_one_sim(dat, "CITS_0.4_constant")
    res_spill3  <- fit_one_sim(dat, "CITS_spill_3pct")
    
    tibble(
      Method    = c("Trd", "CITS", "CITS_0.4_constant", "CITS_spill_3pct"),
      estimate  = c(as.numeric(res_trd$estimate)[1],
                    as.numeric(res_cits$estimate)[1],
                    as.numeric(res_cits04$estimate)[1],
                    as.numeric(res_spill3$estimate)[1]),
      se        = c(as.numeric(res_trd$se)[1],
                    as.numeric(res_cits$se)[1],
                    as.numeric(res_cits04$se)[1],
                    as.numeric(res_spill3$se)[1]),
      pvalue    = c(as.numeric(res_trd$p)[1],
                    as.numeric(res_cits$p)[1],
                    as.numeric(res_cits04$p)[1],
                    as.numeric(res_spill3$p)[1]),
      converged = c(as.logical(res_trd$converged)[1],
                    as.logical(res_cits$converged)[1],
                    as.logical(res_cits04$converged)[1],
                    as.logical(res_spill3$converged)[1]),
      true_B    = NA_real_
    )
  }) %>%
  ungroup() %>%
  mutate(true_B = true_B_fun(as.character(EstimandScenario)))

## Quick check:
dplyr::count(estimates, rho, n, EstimandScenario, Method)

## --- 6. Performance summary across simulations ---
perf <- estimates %>%
  group_by(rho, n, EstimandScenario, Method) %>%
  summarise(
    ## Total sims and convergence info
    R_total       = n(),
    R_converged   = sum(converged),
    nonconv_rate  = ((R_total - R_converged) / R_total) * 100,
    
    ## True value (constant within group; restrict to converged if possible)
    true_val = dplyr::first(true_B[converged %in% TRUE]),
    B_hat    = round(mean(estimate[converged], na.rm = TRUE), 4),
    
    ## ---- Bias & its MCSE ----
    Bias = if (R_converged > 0) {
      round(mean(estimate[converged], na.rm = TRUE) - dplyr::first(true_B[converged]), 4)
    } else NA_real_,
    
    MCSE_Bias = if (R_converged > 1) {
      est  <- estimate[converged]
      mhat <- mean(est)
      round(sqrt(sum((est - mhat)^2) / (R_converged * (R_converged - 1))), 4)
    } else NA_real_,
    
    Bias_percent = if (!is.na(true_val) && true_val != 0 && R_converged > 0) {
      round(100 * abs(Bias / true_val), 2)
    } else NA_real_,
    
    ## ---- MSE & its MCSE ----
    MSE_estimate = if (R_converged > 0) {
      round(mean((estimate[converged] - true_B[converged])^2), 4)
    } else NA_real_,
    
    MCSE_MSE = if (R_converged > 1) {
      sq_err <- (estimate[converged] - true_B[converged])^2
      m      <- mean(sq_err)
      sqrt(sum((sq_err - m)^2) / (R_converged * (R_converged - 1)))
    } else NA_real_,
    
    ## ---- Coverage & MCSE (95% Wald CI) ----
    coverage = if (R_converged > 0) {
      mean(
        (estimate[converged] - 1.96 * se[converged]) <= true_B[converged] &
          true_B[converged] <= (estimate[converged] + 1.96 * se[converged])
      )
    } else NA_real_,
    
    MCSE_coverage = if (R_converged > 0 && !is.na(coverage)) {
      round(sqrt(coverage * (1 - coverage) / R_converged), 4)
    } else NA_real_,
    
    ## ---- Power & MCSE (H0: no effect) ----
    power = if (R_converged > 0) {
      mean(pvalue[converged] <= 0.05, na.rm = TRUE)
    } else NA_real_,
    
    MCSE_power = if (R_converged > 0 && !is.na(power)) {
      round(sqrt(power * (1 - power) / R_converged), 4)
    } else NA_real_,
    
    ## ---- Empirical vs model-based SE ----
    Empirical_SE      = if (R_converged > 1) round(sd(estimate[converged]), 4) else NA_real_,
    avg_model_SE      = if (R_converged > 0) round(sqrt(mean((se[converged])^2)), 4) else NA_real_,
    ratio_Emp_ModelSE = if (!is.na(Empirical_SE) && !is.na(avg_model_SE) && avg_model_SE > 0) {
      round(Empirical_SE / avg_model_SE, 4)
    } else NA_real_,
    
    .groups = "drop"
  ) %>%
  mutate(
    power_pct      = round(100 * power, 2),
    MCSE_power_pct = round(100 * MCSE_power, 4),
    coverage_pct   = round(100 * coverage, 2),
    MCSE_cov_pct   = round(100 * MCSE_coverage, 4)
  )

head(perf)

## --- 7. Save simulation-level estimates and performance summary ---
if (!dir.exists("data")) dir.create("data", recursive = TRUE)

saveRDS(estimates, "data/sim_estimates.rds")
write.csv(estimates, "data/sim_estimates.csv", row.names = FALSE)

saveRDS(perf, "data/performance_summary.rds")
write.csv(perf, "data/performance_summary.csv", row.names = FALSE)
##monitoring run
#tail -f ml_analysis_runtime.log

