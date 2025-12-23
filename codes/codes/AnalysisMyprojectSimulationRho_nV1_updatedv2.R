rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

library(MASS)
library(dplyr)
library(sandwich)
library(parallel)

## =========================
## USER SETTINGS
## =========================
N_CORES_WIN  <- 2   # Windows: PSOCK parallel (duplicates memory)
N_CORES_UNIX <- 4   # Linux/Lengau: fork parallel (shared memory)
SEED         <- 123 # optional reproducibility

delta  <- 1e-4
fam_nb <- MASS::negative.binomial(350)

true_B_fun <- function(scen) {
  dplyr::case_when(
    scen == "Small"    ~ -0.0408,
    scen == "Moderate" ~ -0.3567,
    scen == "Large"    ~ -0.5108,
    TRUE               ~ NA_real_
  )
}

fit_one_sim <- function(dat, method = c("Trd", "CITS", "CITS_0.4_constant")) {
  method <- match.arg(method)
  L <- floor(nrow(dat)^(1/4))
  
  fit <- tryCatch(
    {
      if (method == "CITS") {
        glm(Y ~ P + offset(log(Z + delta)), family = fam_nb, data = dat)
      } else if (method == "CITS_0.4_constant") {
        glm(Y ~ P + offset(log(Z_zeropoint4Cor + delta)), family = fam_nb, data = dat)
      } else {
        glm(Y ~ t + P + X1, family = fam_nb, data = dat)
      }
    },
    error = function(e) NULL
  )
  
  if (is.null(fit) || !isTRUE(fit$converged)) {
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

## =========================
## 1) Load big data
## =========================
AllRhoCombined <- readRDS("data/AllRhoCombined_diffAR_ZvsY.rds")

## =========================
## 2) Keep only needed columns (memory saver)
## =========================
keep_cols <- c("rho","n","EstimandScenario","j","t","t0","P","Y","X1","Z","Z_zeropoint4Cor")
AllRhoCombined <- AllRhoCombined[, keep_cols]
gc()

## =========================
## 3) Order once
## =========================
AllRhoCombined <- AllRhoCombined %>%
  arrange(rho, n, EstimandScenario, j, t)

## =========================
## 4) Build aligned grouping objects (CRITICAL to match non-parallel)
## =========================
gdf <- AllRhoCombined %>%
  group_by(rho, n, EstimandScenario, j) %>%
  mutate(.gid = dplyr::cur_group_id()) %>%
  ungroup()

keys_tbl <- gdf %>%
  distinct(.gid, rho, n, EstimandScenario, j) %>%
  arrange(.gid)

idx_list <- split(seq_len(nrow(gdf)), gdf$.gid)
stopifnot(length(idx_list) == nrow(keys_tbl))

## =========================
## 5) Worker function (models unchanged)
## =========================
worker_fun <- function(i, idx_list, keys_tbl, gdf) {
  ii  <- idx_list[[i]]
  dat <- gdf[ii, , drop = FALSE]
  
  res_trd    <- fit_one_sim(dat, "Trd")
  res_cits   <- fit_one_sim(dat, "CITS")
  res_cits04 <- fit_one_sim(dat, "CITS_0.4_constant")
  
  data.frame(
    rho              = keys_tbl$rho[i],
    n                = keys_tbl$n[i],
    EstimandScenario = keys_tbl$EstimandScenario[i],
    j                = keys_tbl$j[i],
    Method           = c("Trd", "CITS", "CITS_0.4_constant"),
    estimate         = c(res_trd$estimate,  res_cits$estimate,  res_cits04$estimate),
    se               = c(res_trd$se,        res_cits$se,        res_cits04$se),
    pvalue           = c(res_trd$p,         res_cits$p,         res_cits04$p),
    converged        = c(res_trd$converged, res_cits$converged, res_cits04$converged),
    stringsAsFactors = FALSE
  )
}

## =========================
## 6) Parallel execution (Windows + Linux/Lengau)
## =========================
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

if (.Platform$OS.type == "unix") {
  
  N_CORES <- N_CORES_UNIX
  message("UNIX detected. Using fork parallel with ", N_CORES, " cores.")
  
  est_list <- parallel::mclapply(
    X = seq_along(idx_list),
    FUN = worker_fun,
    idx_list = idx_list,
    keys_tbl = keys_tbl,
    gdf = gdf,
    mc.cores = N_CORES,
    mc.preschedule = FALSE
  )
  
} else {
  
  N_CORES <- N_CORES_WIN
  message("Windows detected. Using PSOCK parallel with ", N_CORES,
          " cores (memory duplication is expected).")
  
  cl <- parallel::makeCluster(N_CORES)
  
  parallel::clusterEvalQ(cl, {
    library(MASS)
    library(dplyr)
    library(sandwich)
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "idx_list","keys_tbl","gdf",
      "delta","fam_nb","true_B_fun","fit_one_sim","worker_fun"
    ),
    envir = environment()
  )
  
  est_list <- parallel::parLapply(
    cl,
    X = seq_along(idx_list),
    fun = worker_fun,
    idx_list = idx_list,
    keys_tbl = keys_tbl,
    gdf = gdf
  )
  
  parallel::stopCluster(cl)
}

estimates <- dplyr::bind_rows(est_list) %>%
  dplyr::mutate(true_B = true_B_fun(as.character(EstimandScenario)))

## =========================
## 7) Performance summary (FORMULAS UNCHANGED)
## =========================
perf <- estimates %>%
  group_by(rho, n, EstimandScenario, Method) %>%
  summarise(
    R_total       = n(),
    R_converged   = sum(converged),
    nonconv_rate  = ((R_total - R_converged) / R_total) * 100,
    true_val      = first(true_B[converged %in% TRUE]),
    B_hat         = round(mean(estimate[converged]), 4),
    
    Bias = if (R_converged > 0) round(mean(estimate[converged]) - first(true_B[converged]), 4) else NA_real_,
    
    MCSE_Bias = if (R_converged > 1) {
      est  <- estimate[converged]
      mhat <- mean(est)
      round(sqrt(sum((est - mhat)^2) / (R_converged * (R_converged - 1))), 4)
    } else NA_real_,
    
    Bias_percent = if (!is.na(true_val) && true_val != 0 && R_converged > 0) {
      round(100 * abs(Bias / true_val), 2)
    } else NA_real_,
    
    MSE_estimate = if (R_converged > 0) round(mean((estimate[converged] - true_B[converged])^2), 4) else NA_real_,
    
    MCSE_MSE = if (R_converged > 1) {
      sq_err <- (estimate[converged] - true_B[converged])^2
      m <- mean(sq_err)
      sqrt(sum((sq_err - m)^2) / (R_converged * (R_converged - 1)))
    } else NA_real_,
    
    coverage = if (R_converged > 0) {
      mean((estimate[converged] - 1.96 * se[converged]) <= true_B[converged] &
             true_B[converged] <= (estimate[converged] + 1.96 * se[converged]))
    } else NA_real_,
    
    MCSE_coverage = if (R_converged > 0 && !is.na(coverage)) round(sqrt(coverage * (1 - coverage) / R_converged), 4) else NA_real_,
    
    power = if (R_converged > 0) mean(pvalue[converged] <= 0.05) else NA_real_,
    MCSE_power = if (R_converged > 0 && !is.na(power)) round(sqrt(power * (1 - power) / R_converged), 4) else NA_real_,
    
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

## =========================
## 8) Save outputs
## =========================
if (!dir.exists("data")) dir.create("data", recursive = TRUE)
saveRDS(estimates, "data/sim_estimates.rds")
write.csv(estimates, "data/sim_estimates.csv", row.names = FALSE)
saveRDS(perf, "data/performance_summary.rds")
write.csv(perf, "data/performance_summary.csv", row.names = FALSE)

## =========================
## 9) Row-count diagnostics
## =========================
cat("\n===== ROW COUNTS / DIAGNOSTICS =====\n")
cat("Rows in gdf (All simulations rows)      =", nrow(gdf), "\n")
cat("Unique simulations (unique gid)        =", nrow(keys_tbl), "\n")
cat("Rows in estimates (should be gid*3)    =", nrow(estimates), "\n")
cat("Rows in perf (groups in summary table) =", nrow(perf), "\n")
cat("Check: gid*3                            =", nrow(keys_tbl) * 3, "\n")
cat("=====================================\n")
