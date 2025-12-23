## =========================================================
##  BLOCK 1: Build performance summary from AllRhoCombined
## =========================================================

rm(list = ls(all.names = TRUE), envir = .GlobalEnv)
#rm(list = ls(all.names = TRUE), envir = .GlobalEnv)
## --- Packages ---
#if (!requireNamespace("MASS", quietly = TRUE))      install.packages("MASS")
#if (!requireNamespace("dplyr", quietly = TRUE))     install.packages("dplyr")
#if (!requireNamespace("sandwich", quietly = TRUE))  install.packages("sandwich")
#if (!requireNamespace("lmtest", quietly = TRUE))    install.packages("lmtest")
#if (!requireNamespace("purrr", quietly = TRUE))     install.packages("purrr")
#if (!requireNamespace("ggplot2", quietly = TRUE))   install.packages("ggplot2")

library(MASS)
library(dplyr)
library(sandwich)
library(lmtest)
library(purrr)
library(ggplot2)

## --- (Optional) Source your own bias/power helpers here ---
## source("R/perf_functions.R")  # e.g. bias(), power_fun(), coverage_fun(), etc.

## --- 1. Load your big simulated data set ---

## OPTION A: from RDS (if you saved it like before) AllRhoCombined_diffAR_ZvsY
#AllRhoCombined <- readRDS("data/AllRhoCombined.rds")
AllRhoCombined <- readRDS("data/AllRhoCombined_diffAR_ZvsY.rds")
## OPTION B: from CSV (uncomment if you’re using a .csv)
# AllRhoCombined <- read.csv("data/AllRhoCombined.csv")

str(AllRhoCombined)

## --- 2. Make sure EstimandScenario is a labelled factor ---

#AllRhoCombined <- AllRhoCombined %>%
  #mutate(
    #EstimandScenario = factor(
      #EstimandScenario,
      #levels = c(1, 2, 3),
      #abels = c("Expected", "minimum", "Maximum")
    #)
  #)

## --- 3. True policy effects for each scenario (on log scale) ---
## These are the B2 values you gave:
## small  = -0.0408
## medium = -0.3567
## large  = -0.5108

true_B_fun <- function(scen) {
  dplyr::case_when(
    scen == "Small"  ~ -0.0408,  # small
    scen == "Moderate" ~ -0.3567,  # moderate
    scen == "Large"  ~ -0.5108,  # large
    TRUE               ~ NA_real_
  )
}

## --- 4. Function to fit one simulation for a given method ---

delta  <- 1e-4                          # avoid log(0)
fam_nb <- MASS::negative.binomial(350)  # same dispersion as simulation

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  as.numeric(cor(x[ok], y[ok]))
}


common_trend_p <- function(dat, zvar = c("Z", "Z_zeropoint4Cor"),
                           delta = 1e-4,
                           fam_nb = MASS::negative.binomial(350)) {
  
  zvar <- match.arg(zvar)
  
  datPre <- dat %>% dplyr::filter(t < t0)
  if (nrow(datPre) < 5) return(NA_real_)
  
  form <- if (zvar == "Z") {
    Y ~ t + offset(log(Z + delta))
  } else {
    Y ~ t + offset(log(Z_zeropoint4Cor + delta))
  }
  
  fitPre <- tryCatch(glm(form, family = fam_nb, data = datPre),
                     error = function(e) NULL)
  
  if (is.null(fitPre) || !isTRUE(fitPre$converged)) return(NA_real_)
  
  # Newey–West HAC (same lag rule as your main code)
  L <- floor(nrow(datPre)^(1/4))
  
  out <- tryCatch({
    V  <- sandwich::NeweyWest(fitPre, lag = L, prewhite = FALSE, adjust = TRUE)
    est <- coef(fitPre)["t"]
    se  <- sqrt(diag(V)["t"])
    z   <- est / se
    2 * pnorm(-abs(z))
  }, error = function(e) NA_real_)
  
  as.numeric(out)
}




fit_one_sim <- function(dat, method = c("Trd", "CITS", "CITS_0.4_constant")) {
  method <- match.arg(method)
  L <- floor(nrow(dat)^(1/4))   # HAC lag as in your code
  
  ## Fit model with basic error protection, track convergence
  fit <- tryCatch(
    {
      if (method == "CITS") {
        ## CITS: Y ~ P + offset(log(Z))
        glm(Y ~ P + offset(log(Z + delta)),
            family = fam_nb, data = dat)
        
      } else if (method == "CITS_0.4_constant") {
        ## CITS with constant AR(1) Z (rho=0.4)
        glm(Y ~ P + offset(log(Z_zeropoint4Cor + delta)),
            family = fam_nb, data = dat)
        
      } else {
        ## Traditional regression: Y ~ t + P + X1
        glm(Y ~ t + P + X1,
            family = fam_nb, data = dat)
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    ## model failed to fit
    return(list(
      estimate  = NA_real_,
      se        = NA_real_,
      p         = NA_real_,
      converged = FALSE
    ))
  }
  
  conv_flag <- isTRUE(fit$converged)
  
  if (!conv_flag) {
    ## no HAC etc. if not converged
    return(list(
      estimate  = NA_real_,
      se        = NA_real_,
      p         = NA_real_,
      converged = FALSE
    ))
  }
  
  ## Newey–West HAC variance and inference, safely
  out <- tryCatch(
    {
      V <- sandwich::NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE)
      
      est <- coef(fit)["P"]
      se  <- sqrt(diag(V)["P"])
      z   <- est / se
      p   <- 2 * pnorm(-abs(z))
      
      list(
        estimate  = est,
        se        = se,
        p         = p,
        converged = TRUE
      )
    },
    error = function(e) {
      ## glm estimated P, but HAC variance failed
      est <- coef(fit)["P"]
      list(
        estimate  = est,        # keep the point estimate
        se        = NA_real_,   # no valid HAC SE
        p         = NA_real_,   # so no valid p-value
        converged = FALSE       # mark as non-converged for perf summaries
      )
    }
  )
  
  out
}

names(AllRhoCombined)
## --- 5. Get simulation–level results for every (rho, n, scenario, j, method) ---

## --- 5. Get simulation–level results for every (rho, n, scenario, j, method) ---
estimates <- AllRhoCombined %>%
  arrange(rho, n, EstimandScenario, j, t) %>%
  group_by(rho, n, EstimandScenario, j) %>%
  group_modify(~{
    dat   <- .x
    
    res_trd   <- fit_one_sim(dat, "Trd")
    res_cits  <- fit_one_sim(dat, "CITS")
    res_cits04 <- fit_one_sim(dat, "CITS_0.4_constant")
    p_valueZvar      <- common_trend_p(dat, zvar = "Z",              delta = delta, fam_nb = fam_nb)
    p_valueZconstant <- common_trend_p(dat, zvar = "Z_zeropoint4Cor",delta = delta, fam_nb = fam_nb)
    cor_YZ <- safe_cor(dat$Y, dat$Z)
    
    # innovations correlation: only if present in AllRhoCombined
    cor_eY_eZ <- if (all(c("e_y","e_z") %in% names(dat))) {
      safe_cor(dat$e_y, dat$e_z)
    } else {
      NA_real_
    }
    
    tibble(
      Method    = c("Trd", "CITS", "CITS_0.4_constant"),
      estimate  = c(as.numeric(res_trd$estimate)[1],
                    as.numeric(res_cits$estimate)[1],
                    as.numeric(res_cits04$estimate)[1]),
      se        = c(as.numeric(res_trd$se)[1],
                    as.numeric(res_cits$se)[1],
                    as.numeric(res_cits04$se)[1]),
      pvalue    = c(as.numeric(res_trd$p)[1],
                    as.numeric(res_cits$p)[1],
                    as.numeric(res_cits04$p)[1]),
      converged = c(as.logical(res_trd$converged)[1],
                    as.logical(res_cits$converged)[1],
                    as.logical(res_cits04$converged)[1]),
      true_B    = NA_real_,
      # NEW: pre-trend/common trend check p-values
      p_valueZvar      = c(NA_real_, p_valueZvar,      NA_real_),
      p_valueZconstant = c(NA_real_, NA_real_, p_valueZconstant),
      
      # NEW: correlations (same for the 3 rows within the same sim)
      cor_YZ    = cor_YZ,
      cor_eY_eZ = cor_eY_eZ,
      
      # NEW: indicators for threshold
      cor_YZ_ge_066    = as.integer(!is.na(cor_YZ)    & cor_YZ    >= 0.66),
      cor_eY_eZ_ge_066 = as.integer(!is.na(cor_eY_eZ) & cor_eY_eZ >= 0.66)
    )
  }) %>%
  ungroup() %>%
  mutate(
    true_B = true_B_fun(as.character(EstimandScenario))
  )

#View(estimates)
## Quick check:
dplyr::count(estimates, rho, n, EstimandScenario, Method)

## --- 6. Summarise across simulations (performance table) ---*2430000

perf <- estimates %>%
  group_by(rho, n, EstimandScenario, Method) %>%
  summarise(
    ## Total sims and convergence info
    R_total       = n(),
    R_converged   = sum(converged),
    nonconv_rate  = ((R_total - R_converged) / R_total)*100,
    mean_cor_YZ = round(mean(cor_YZ, na.rm = TRUE), 4),
    
    prop_cor_YZ_ge_066 = round(mean(cor_YZ_ge_066, na.rm = TRUE), 4),
    
    # innovations (will be NA if innovations weren't saved)
    mean_cor_eY_eZ = if (all(is.na(cor_eY_eZ))) NA_real_ else round(mean(cor_eY_eZ, na.rm = TRUE), 4),
    
    prop_cor_eY_eZ_ge_066 = if (all(is.na(cor_eY_eZ_ge_066))) NA_real_ else round(mean(cor_eY_eZ_ge_066, na.rm = TRUE), 4),
    
    ## True value (constant within group, restrict to converged if possible)
    true_val = dplyr::first(true_B[converged %in% TRUE]),
    B_hat=round(mean(estimate[converged]),4),
    
    ## ---- Bias & its MCSE ----
    ## If you prefer your own bias() function, you can replace this block by:
    ##   bias_res <- bias( ... )
    ##   Bias          = bias_res$Bias
    ##   Bias_percent  = bias_res$Bias_percent
    ##
    Bias = if (R_converged > 0) {
      #mean(estimate[converged] - true_B[converged])
      round(mean(estimate[converged]) - dplyr::first(true_B[converged]),4)
    } else {
      NA_real_
    },
    MCSE_Bias = if (R_converged > 1) {
      est  <- estimate[converged]          # θ̂_i
      mhat <- mean(est)                    # \bar{θ̂}
      round(sqrt( sum((est - mhat)^2) / (R_converged * (R_converged - 1)) ),4)
    } else {
      NA_real_
    }
    ,
    Bias_percent = if (!is.na(true_val) && true_val != 0 && R_converged > 0) {
      round(100 * abs(Bias / true_val),2)
    } else {
      NA_real_
    },
    #MCSE_Bias_percent = if (!is.na(true_val) && true_val != 0 && R_converged > 1) {
      #rel_bias_i <- (estimate[converged] - true_B[converged]) / true_val
      #100 * sd(rel_bias_i) / sqrt(R_converged)
    #} else {
      #NA_real_
    #},
    
    ## ---- MSE & its MCSE ----
    MSE_estimate = if (R_converged > 0) {
      round(mean((estimate[converged] - true_B[converged])^2),4)
    } else {
      NA_real_
    },
    MCSE_MSE = if (R_converged > 1) {
      sq_err <- (estimate[converged] - true_B[converged])^2  # X_i
      m      <- mean(sq_err)                                 # MSE = mean(X_i)
      sqrt( sum((sq_err - m)^2) / (R_converged * (R_converged - 1)) )
    } else {
      NA_real_
    },
    
    ## ---- Coverage & MCSE (95% Wald CI) ----
    coverage = if (R_converged > 0) {
      mean(
        (estimate[converged] - 1.96 * se[converged]) <= true_B[converged] &
          true_B[converged] <= (estimate[converged] + 1.96 * se[converged])
      )
    } else {
      NA_real_
    },
    MCSE_coverage = if (R_converged > 0 && !is.na(coverage)) {
      round(sqrt(coverage * (1 - coverage) / R_converged),4)
    } else {
      NA_real_
    },
    
    ## ---- Power & MCSE (H0: no effect) ----
    ## If you have your own power() function, you can replace this block.
    power = if (R_converged > 0) {
      mean(pvalue[converged] <= 0.05)
    } else {
      NA_real_
    },
    MCSE_power = if (R_converged > 0 && !is.na(power)) {
      round(sqrt(power * (1 - power) / R_converged),4)
    } else {
      NA_real_
    },
    pretrend_pass = if (R_converged > 0) { mean( dplyr::case_when( Method == "CITS" ~ (p_valueZvar > 0.05), Method == "CITS_0.4_constant" ~ (p_valueZconstant > 0.05), TRUE ~ NA ), na.rm = TRUE ) } else NA_real_
    ,
    
    ## ---- Empirical vs model-based SE ----
    Empirical_SE      = if (R_converged > 1) round(sd(estimate[converged]),4) else NA_real_,
    avg_model_SE      = if (R_converged > 0) round(sqrt(mean((se[converged])^2)),4)    else NA_real_,
    ratio_Emp_ModelSE = if (!is.na(Empirical_SE) && !is.na(avg_model_SE) && avg_model_SE > 0) {
      round(Empirical_SE / avg_model_SE,4)
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    power_pct        = round(100 * power,2),
    MCSE_power_pct   = round(100 * MCSE_power,4),
    coverage_pct     = round(100 * coverage,2),
    MCSE_cov_pct     = round(100 * MCSE_coverage,4)
  )

## Look at the first rows
head(perf)

## --- 7. Save simulation-level estimates and performance summary ---

if (!dir.exists("data")) dir.create("data", recursive = TRUE)

## Store estimates and their SEs (per-sim, per-method)
saveRDS(estimates, "data/sim_estimates.rds")
write.csv(estimates, "data/sim_estimates.csv", row.names = FALSE)

## Store performance summary
saveRDS(perf, "data/performance_summary.rds")
write.csv(perf, "data/performance_summary.csv", row.names = FALSE)

