
#packages
library(dplyr)
library(patchwork)
#library(ggplot2)
library(grid)
library(ggplot2)
library(tidyr)
if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}
library(ggridges)
##reading data
perf<-readRDS("data/performance_summary.rds")
str(perf)
##creating a table of performance measures to be used in R_shinydashbord
perfRshiny <- perf %>%
  mutate(
    `Policy effect (Empirical SE)`        = sprintf("%.4f (%.4f)", B_hat, Empirical_SE),
    `Bias (MCSE)`                         = sprintf("%.4f (%.4f)", Bias, MCSE_Bias),
    `Coverage 95% (MCSE)`                 = sprintf("%.2f (%.4f)", coverage_pct, MCSE_cov_pct),
    `Power (MCSE)`                        = sprintf("%.2f (%.4f)", power_pct, MCSE_power_pct),
    `Empirical SE (Model-based SE)`       = sprintf("%.4f (%.4f)", Empirical_SE, avg_model_SE),
    `Mean square error (MCSE)`            = sprintf("%.4f (%.4f)", MSE_estimate, MCSE_MSE)
  ) %>%
  select(
    rho, n, EstimandScenario, Method,
    `Policy effect (Empirical SE)`,
    `Bias (MCSE)`,
    `Coverage 95% (MCSE)`,
    `Power (MCSE)`,
    `Empirical SE (Model-based SE)`,
    `Mean square error (MCSE)`,
    true_val
  )


write.csv(perfRshiny,
          "data/perfRshiny.csv",
          row.names = FALSE)
#levels = c(1, 2, 3),
#labels = c("Expected", "minimum", "Maximum"))

perf_plot <- perf %>%
  mutate(
    EstimandScenario = if (!is.factor(EstimandScenario)) {
      factor(EstimandScenario,
             levels = c("Small", "Moderate", "Large"))
    } else {
      EstimandScenario
    },
    
    rho_f = if (!is.factor(rho)) {
      factor(
        rho,
        levels = c("0", "0.2", "0.4", "0.6", "0.8"),
        labels = c("ρ = 0.0", "ρ = 0.2", "ρ = 0.4", "ρ = 0.6", "ρ = 0.8")
      )
    } else {
      rho
    },
    
    Method = if (!is.factor(Method)) {
      factor(Method, levels = c("Trd", "CITS","CITS_0.4_constant"))
    } else {
      Method
    },
    
    ## --- dynamic facet labels using true_val ---
    facet_label = case_when(
      EstimandScenario == "Small"    ~ paste0("Small(", round(true_val, 4), ")"),
      EstimandScenario == "Moderate" ~ paste0("Moderate(", round(true_val, 4), ")"),
      EstimandScenario == "Large"    ~ paste0("Large(", round(true_val, 4), ")"),
      TRUE ~ as.character(EstimandScenario)
    )
  ) %>%
  
  ## ---- force the ordering of facet_label ----
mutate(
  facet_label = factor(
    facet_label,
    levels = c(
      facet_label[EstimandScenario == "Small"][1],
      facet_label[EstimandScenario == "Moderate"][1],
      facet_label[EstimandScenario == "Large"][1]
    )
  )
)


 #STANDARD ERROR PLOTS
##traditional regression plot
p_trd2 <- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(aes(x = n, group = 1)) +
  
  # ---- Empirical SE ----
geom_line(aes(y = Empirical_SE, colour = "Empirical SE")) +
  geom_point(aes(y = Empirical_SE, colour = "Empirical SE"), size = 1.5) +
  
  # ---- Model-based SE ----
geom_line(aes(y = avg_model_SE, colour = "Model-based SE"), linetype = "dashed") +
  geom_point(aes(y = avg_model_SE, colour = "Model-based SE"), size = 1.5) +
  
  facet_grid(facet_label ~ rho_f) +
  scale_colour_manual(
    name = "SE type",
    values = c(
      "Empirical SE" = "#F8766D",   # red (your original Trd colour)
      "Model-based SE" = "#00BFC4" # turquoise to match CITS
    )
  ) +
  labs(
    x = "Length of time series (n)",
    y = "Standard Error(SE)",
    title = "A"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

p_trd2

ggsave(
  filename = "graphs/p_trd_SE_traditional.png",
  plot     = p_trd2,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)

ggsave(
  filename = "graphs/p_trd_SE_traditional.pdf",
  plot     = p_trd2,
  width    = 8,      # inches
  height   = 4.5,
  device   = cairo_pdf
)


##CITS Controlled Interrupted Time Series (CITS)
p_cits2 <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(aes(x = n, group = 1)) +
  
  # ---- Empirical SE (your existing turquoise) ----
geom_line(aes(y = Empirical_SE, colour = "Empirical SE")) +
  geom_point(aes(y = Empirical_SE, colour = "Empirical SE"), size = 1.5) +
  
  # ---- Model-based SE (add second colour, dashed) ----
geom_line(aes(y = avg_model_SE, colour = "Model-based SE"), 
          linetype = "dashed", size = 0.9) +
  geom_point(aes(y = avg_model_SE, colour = "Model-based SE"), size = 1.3) +
  
  facet_grid(facet_label~ rho_f) +
  scale_colour_manual(
    name = "SE Type",
    values = c(
      "Empirical SE"   = "#F8766D",  # your CITS turquoise
      "Model-based SE" = "#00BFC4"  # soft red (contrasts nicely)
    )
  ) +
  labs(
    x = "Length of time series (n)",
    y = "Standard Error(SE)",
    title = "B"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

p_cits2
ggsave(
  filename = "graphs/p_cits_SE.png",
  plot     = p_cits2,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)

ggsave(
  filename = "graphs/p_cits_SE.pdf",
  plot     = p_cits2,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


#original scale
#log scale
perf_plotcitstrd<-perf_plot[perf_plot$Method!="CITS_0.4_constant",]
p_combined_log2 <- ggplot(
  perf_plotcitstrd,
  aes(x = n, y = Empirical_SE, colour = Method, group = Method)
) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(facet_label ~ rho_f) +
  scale_y_log10() +
  labs(
    x = "Length of time series (n)",
    y = "Empirical SE(log scale)",
    colour = "Method"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

p_combined_log2


#p_cits2
ggsave(
  filename = "graphs/p_trd_SE_cits_trdCombined.png",
  plot     = p_combined_log2,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)


ggsave(
  filename = "graphs/p_trd_SE_cits_trdCombined.pdf",
  plot     = p_combined_log2,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


###Power plot Multivariable regression
##traditional regression plot
#rho_f<-c(0,0.2,0.4,0.6,0.8)
p_trd_power <- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(
    aes(
      x      = n,
      y      = power_pct,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Power",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_trd_power

ggsave(
  filename = "graphs/p_trd_Power.png",
  plot     = p_trd_power,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)




#CITS
p_cits_power <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      = power_pct,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Power",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_cits_power

ggsave(
  filename = "graphs/p_cits_power.png",
  plot     = p_cits_power,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)






##Coverage graphs
#cits
p_trd_cits <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      = coverage_pct,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +geom_hline(yintercept = 95, colour = "black", linetype = "solid", size = 0.8)+
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Coverage(%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_trd_cits 

ggsave(
  filename = "graphs/p_cits_coverage.png",
  plot     =p_trd_cits ,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)






#coveragetrd
p_trd_cova<- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(
    aes(
      x      = n,
      y      = coverage_pct,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +geom_hline(yintercept = 95, colour = "black", linetype = "solid", size = 0.8)+
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Coverage(%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_trd_cova 

ggsave(
  filename = "graphs/p_trd_coverage.png",
  plot     =p_trd_cova ,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)











#MSE
##trd MSE
p_trd_mse <- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(
    aes(
      x      = n,
      y      =MSE_estimate,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Mean Squared Error(MSE)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_trd_mse

ggsave(
  filename = "graphs/p_trd_mse.png",
  plot     =p_trd_mse,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)


#cits
p_cits_mse <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      =MSE_estimate,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Mean Squared Error(MSE)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_cits_mse

ggsave(
  filename = "graphs/p_cits_mse.png",
  plot     =p_cits_mse,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)







##Plotting point estimates
## If you saved estimates earlier:
## 1. Base data (what you already have)
estimates <- readRDS("data/sim_estimates.rds")
##Check som summaries
mean(estimates$estimate[estimates$rho==0.8 & estimates$n==100 & estimates$EstimandScenario=="Large"])




est_plot <- estimates %>%
  filter(converged, is.finite(estimate)) %>%
  mutate(
    EstimandScenario = if (!is.factor(EstimandScenario)) {
      factor(EstimandScenario,
             levels = c("Small", "Moderate", "Large"))
    } else EstimandScenario,
    
    rho_char = as.character(rho),
    rho_f = factor(
      rho_char,
      levels = c("0", "0.2", "0.4", "0.6", "0.8"),
      labels = c("ρ = 0.0", "ρ = 0.2", "ρ = 0.4", "ρ = 0.6", "ρ = 0.8")
    ),
    n = if (!is.numeric(n)) as.numeric(as.character(n)) else n
  ) %>%
  dplyr::select(-rho_char)

## 2. Scenario labels + true effects (fix facet order here)
scenario_labels <- est_plot %>%
  distinct(EstimandScenario, true_B) %>%
  arrange(EstimandScenario) %>%   # ensures Small, Moderate, Large
  mutate(
    facet_label = paste0(EstimandScenario, " (", round(true_B, 4), ")"),
    facet_label = factor(facet_label, levels = facet_label)
  )

## This will also be used for the dashed lines
true_lines <- scenario_labels

## 3. Subset to selected n and join facet labels
selected_n <- c(12, 24, 40, 60, 80, 100)

est_trd <- est_plot %>%
  filter(Method == "Trd", n %in% selected_n) %>%
  left_join(
    scenario_labels %>% select(EstimandScenario, facet_label),
    by = "EstimandScenario"
  ) %>%
  mutate(
    n_f = factor(n, levels = selected_n)   # 12 at bottom, 100 at top
  )

## 4. Plot
p_trd_ridge <- ggplot(
  est_trd,
  aes(
    x      = estimate,
    y      = n_f,
    colour = rho_f,
    group  = interaction(n_f, rho_f)
  )
) +
  geom_density_ridges(
    scale          = 1.2,
    rel_min_height = 0.001,
    size           =1,
    alpha          = 0.7
  ) +
  geom_vline(
    data        = true_lines,
    aes(xintercept = true_B),
    linetype    = "dashed",
    colour      = "black",
    inherit.aes = FALSE
  ) +
  facet_grid(. ~ facet_label) +
  labs(
    title  = "A",
    x      = "Estimated intervention effect (log IRR)",
    y      = "Length of time series (n)",
    colour = "Autocorrelation (ρ)"
  ) +
  ggridges::theme_ridges() +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(size = 11, face = "bold"),  # smaller labels
    panel.grid.minor = element_blank(),
    axis.title.x     = element_text(hjust = 0.5, vjust = -0.5) 
  )
# optional: tighten x-range if you want
# + coord_cartesian(xlim = c(-2, 1))

p_trd_ridge

ggsave(
  "graphs/Trd_ridgeline_estimates.png",
  p_trd_ridge,
  width  = 8,
  height = 4.5,
  dpi    = 300
)

ggsave(
  filename = "graphs/Trd_ridgeline_estimates.pdf",
  plot     = p_trd_ridge,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)

##CITS
est_plot <- estimates %>%
  filter(converged, is.finite(estimate)) %>%
  mutate(
    EstimandScenario = if (!is.factor(EstimandScenario)) {
      factor(EstimandScenario,
             levels = c("Small", "Moderate", "Large"))
    } else EstimandScenario,
    
    rho_char = as.character(rho),
    rho_f = factor(
      rho_char,
      levels = c("0", "0.2", "0.4", "0.6", "0.8"),
      labels = c("ρ = 0.0", "ρ = 0.2", "ρ = 0.4", "ρ = 0.6", "ρ = 0.8")
    ),
    n = if (!is.numeric(n)) as.numeric(as.character(n)) else n
  ) %>%
  select(-rho_char)

## 2. Scenario labels + true effects (fix facet order here)
scenario_labels <- est_plot %>%
  distinct(EstimandScenario, true_B) %>%
  arrange(EstimandScenario) %>%   # ensures Small, Moderate, Large
  mutate(
    facet_label = paste0(EstimandScenario, " (", round(true_B, 4), ")"),
    facet_label = factor(facet_label, levels = facet_label)
  )

## This will also be used for the dashed lines
true_lines <- scenario_labels

## 3. Subset to selected n and join facet labels
selected_n <- c(12, 24, 40, 60, 80, 100)

est_trd <- est_plot %>%
  filter(Method == "CITS", n %in% selected_n) %>%
  left_join(
    scenario_labels %>% select(EstimandScenario, facet_label),
    by = "EstimandScenario"
  ) %>%
  mutate(
    n_f = factor(n, levels = selected_n)   # 12 at bottom, 100 at top
  )

## 4. Plot
p_trd_ridge <- ggplot(
  est_trd,
  aes(
    x      = estimate,
    y      = n_f,
    colour = rho_f,
    group  = interaction(n_f, rho_f)
  )
) +
  geom_density_ridges(
    scale          = 1.2,
    rel_min_height = 0.001,
    size           =1,
    alpha          = 0.7
  ) +
  geom_vline(
    data        = true_lines,
    aes(xintercept = true_B),
    linetype    = "dashed",
    colour      = "black",
    inherit.aes = FALSE
  ) +
  facet_grid(. ~ facet_label) +
  labs(
    title  = "B",
    x      = "Estimated intervention effect (log IRR)",
    y      = "Length of time series (n)",
    colour = "Autocorrelation (ρ)"
  ) +
  ggridges::theme_ridges() +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(size = 11, face = "bold"),  # smaller labels
    panel.grid.minor = element_blank(),
    axis.title.x     = element_text(hjust = 0.5, vjust = -0.5)
  )
# optional: tighten x-range if you want
# + coord_cartesian(xlim = c(-2, 1))

p_trd_ridge

ggsave(
  "graphs/CITS_ridgeline_estimates.png",
  p_trd_ridge,
  width  = 8,
  height = 4.5,
  dpi    = 300
)

ggsave(
  "graphs/CITS_ridgeline_estimates.pdf",
  p_trd_ridge,
  width  = 8,
  height = 4.5,
  dpi    = 300,
  device   = cairo_pdf
)




##Bias
p_trd_bias <- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(
    aes(
      x      = n,
      y      = Bias,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Bias",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "A"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_trd_bias 

ggsave(
  filename = "graphs/p_trd_bias.png",
  plot     = p_trd_bias ,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)




#cits
p_cits_bias <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      = Bias,
      colour = rho_f,     # colour by rho (nice formatted labels)
      shape  = rho_f,     # different point shapes per rho
      group  = rho_f      # one line per rho
    )
  ) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(facet_label ~ .) +   # 3 rows, 1 column
  labs(
    x = "Length of time series (n)",
    y = "Bias",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "B"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p_cits_bias 

ggsave(
  filename = "graphs/p_cits_bias.png",
  plot     = p_cits_bias ,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)


##lets create some table 
ns_to_show<-c(c(12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,80,100,150))
perf_clean <- perf %>%
  mutate(
    EstimandScenario = if (!is.factor(EstimandScenario)) {
      factor(EstimandScenario,
             levels = c("Small", "Moderate", "Large"))
    } else EstimandScenario,
    
    rho = if (!is.factor(rho)) {
      factor(
        rho,
        levels = c(0, 0.2, 0.4, 0.6, 0.8),
        labels = c("0.0", "0.2", "0.4", "0.6", "0.8")
      )
    } else rho,
    
    Method = if (!is.factor(Method)) {
      factor(Method, levels = c("Trd", "CITS","CITS_0.4_constant"))
    } else Method,
    
    # dynamic labels with true_val (same idea as your facet_label)
    effect_label = case_when(
      EstimandScenario == "Small"    ~ paste0("Small(",    round(true_val, 4), ")"),
      EstimandScenario == "Moderate" ~ paste0("Moderate(", round(true_val, 4), ")"),
      EstimandScenario == "Large"    ~ paste0("Large(",    round(true_val, 4), ")"),
      TRUE ~ as.character(EstimandScenario)
    )
  )

table_bias_empse <- perf_clean %>%
  filter(n %in% ns_to_show) %>%# ---- FIX METHOD LABEL HERE ----
mutate(
  Method = case_when(
    Method == "CITS"              ~ "CITS",
    Method == "CITS_0.4_constant" ~ "CITS(rho = 0.4)",
    Method == "Trd"               ~ "Trd",
    TRUE                          ~ as.character(Method)
  )
) %>%
  
  # string for each cell: Estimate (Empirical_SE)
  mutate(
    cell = sprintf("%.4f (%.4f)", B_hat, Empirical_SE),
    
    # column key = Effect-with-true-val + rho
    col_key = paste0(effect_label, "_rho", rho)
  ) %>%
  
  select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small(-0.0408)", "Moderate(-0.3567)", "Large(-0.5108)")

col_order <- c(
  "n", "Method",
  unlist(lapply(effects, function(eff) paste0(eff, "_rho", rho_seq)))
)

estimates <- table_bias_empse[, col_order]

# Row 1: effect names (repeated for each rho)
`Effect size`<- c(
  "Length of series (n)",          # col 1
  "Method",                        # col 2
  rep("Small effect (-0.0408)",    length(rho_seq)),
  rep("Moderate effect (-0.3567)", length(rho_seq)),
  rep("Large effect (-0.5108)",    length(rho_seq))
)

# Row 2: rho values under each effect
`Correlation`<- c(
  "",                              # under n
  "",                              # under Method
  rep(rho_seq, times = length(effects))
)

# Coerce everything to character so we can rbind
estimates_chr <- estimates %>% mutate(across(everything(), as.character))

pretty_table <- rbind(
  `Effect size`,
  `Correlation`,
  as.matrix(estimates_chr)
)

# Optional: remove column names so Excel doesn’t show them
colnames(pretty_table) <- rep("", ncol(pretty_table))
View(pretty_table)
write.csv(pretty_table,
          "data/table_B_hatEstimate_empse_pretty.csv",
          row.names = FALSE)
write.csv(pretty_table,
          "data/point_estimate.csv",
          row.names = FALSE)


##Thanks nice table

##Bias
table_bias_empse <- perf_clean %>%
  filter(n %in% ns_to_show) %>%# ---- FIX METHOD LABEL HERE ----
mutate(
  Method = case_when(
    Method == "CITS"              ~ "CITS",
    Method == "CITS_0.4_constant" ~ "CITS(rho = 0.4)",
    Method == "Trd"               ~ "Trd",
    TRUE                          ~ as.character(Method)
  )
) %>%
  
  # string for each cell: Estimate (Empirical_SE)
  mutate(
    cell = sprintf("%.4f (%.4f)", Bias,MCSE_Bias),
    
    # column key = Effect-with-true-val + rho
    col_key = paste0(effect_label, "_rho", rho)
  ) %>%
  
  select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small(-0.0408)", "Moderate(-0.3567)", "Large(-0.5108)")

col_order <- c(
  "n", "Method",
  unlist(lapply(effects, function(eff) paste0(eff, "_rho", rho_seq)))
)

estimates <- table_bias_empse[, col_order]

# Row 1: effect names (repeated for each rho)
`Effect size`<- c(
  "Length of series (n)",          # col 1
  "Method",                        # col 2
  rep("Small effect (-0.0408)",    length(rho_seq)),
  rep("Moderate effect (-0.3567)", length(rho_seq)),
  rep("Large effect (-0.5108)",    length(rho_seq))
)

# Row 2: rho values under each effect
`Correlation`<- c(
  "",                              # under n
  "",                              # under Method
  rep(rho_seq, times = length(effects))
)

# Coerce everything to character so we can rbind
estimates_chr <- estimates %>% mutate(across(everything(), as.character))

pretty_table <- rbind(
  `Effect size`,
  `Correlation`,
  as.matrix(estimates_chr)
)

# Optional: remove column names so Excel doesn’t show them
colnames(pretty_table) <- rep("", ncol(pretty_table))

write.csv(pretty_table,
          "data/table_B_hatEstimate_empse_pretty.csv",
          row.names = FALSE)
write.csv(pretty_table,
          "data/bias.csv",
          row.names = FALSE)




##Coverage table
table_bias_empse <- perf_clean %>%
  filter(n %in% ns_to_show) %>%
  
  # ---- FIX METHOD LABEL HERE ----
mutate(
  Method = case_when(
    Method == "CITS"              ~ "CITS",
    Method == "CITS_0.4_constant" ~ "CITS(rho = 0.4)",
    Method == "Trd"               ~ "Trd",
    TRUE                          ~ as.character(Method)
  )
) %>%
  
  # string for each cell: Coverage (MCSE)
  mutate(
    cell = sprintf("%.2f (%.4f)", coverage_pct, MCSE_cov_pct),
    col_key = paste0(effect_label, "_rho", rho)
  ) %>%
  
  select(n, Method, col_key, cell) %>%
  distinct() %>%
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small(-0.0408)", "Moderate(-0.3567)", "Large(-0.5108)")

col_order <- c(
  "n", "Method",
  unlist(lapply(effects, function(eff) paste0(eff, "_rho", rho_seq)))
)

estimates <- table_bias_empse[, col_order]

# Row 1: effect names (repeated for each rho)
`Effect size`<- c(
  "Length of series (n)",          # col 1
  "Method",                        # col 2
  rep("Small effect (-0.0408)",    length(rho_seq)),
  rep("Moderate effect (-0.3567)", length(rho_seq)),
  rep("Large effect (-0.5108)",    length(rho_seq))
)

# Row 2: rho values under each effect
`Correlation`<- c(
  "",                              # under n
  "",                              # under Method
  rep(rho_seq, times = length(effects))
)

# Coerce everything to character so we can rbind
estimates_chr <- estimates %>% mutate(across(everything(), as.character))

pretty_table <- rbind(
  `Effect size`,
  `Correlation`,
  as.matrix(estimates_chr)
)

# Optional: remove column names so Excel doesn’t show them
colnames(pretty_table) <- rep("", ncol(pretty_table))
View(pretty_table)
write.csv(pretty_table,
          "data/table_COV_pct_mcse_pretty.csv",
          row.names = FALSE)

#table(pretty_table)

write.csv(pretty_table,
          "data/coverage.csv",
          row.names = FALSE)

##Thanks nice table 



##Power
table_bias_empse <- perf_clean %>%
  filter(n %in% ns_to_show) %>% # ---- FIX METHOD LABEL HERE ----
mutate(
  Method = case_when(
    Method == "CITS"              ~ "CITS",
    Method == "CITS_0.4_constant" ~ "CITS(rho = 0.4)",
    Method == "Trd"               ~ "Trd",
    TRUE                          ~ as.character(Method)
  )
)%>%
  
  # string for each cell: Estimate (Empirical_SE)
  mutate(
    cell = sprintf("%.2f (%.4f)", power_pct, MCSE_power_pct),
    
    # column key = Effect-with-true-val + rho
    col_key = paste0(effect_label, "_rho", rho)
  ) %>%
  
  select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small(-0.0408)", "Moderate(-0.3567)", "Large(-0.5108)")

col_order <- c(
  "n", "Method",
  unlist(lapply(effects, function(eff) paste0(eff, "_rho", rho_seq)))
)

estimates <- table_bias_empse[, col_order]

# Row 1: effect names (repeated for each rho)
 `Effect size`<- c(
  "Length of series (n)",          # col 1
  "Method",                        # col 2
  rep("Small effect (-0.0408)",    length(rho_seq)),
  rep("Moderate effect (-0.3567)", length(rho_seq)),
  rep("Large effect (-0.5108)",    length(rho_seq))
)

# Row 2: rho values under each effect
`Correlation`<- c(
  "",                              # under n
  "",                              # under Method
  rep(rho_seq, times = length(effects))
)

# Coerce everything to character so we can rbind
estimates_chr <- estimates %>% mutate(across(everything(), as.character))

pretty_table <- rbind(
  `Effect size`,
  `Correlation`,
  as.matrix(estimates_chr)
)

# Optional: remove column names so Excel doesn’t show them
colnames(pretty_table) <- rep("", ncol(pretty_table))

View(pretty_table)
write.csv(pretty_table,
          "data/table_POWER_pct_mcse_pretty.csv",
          row.names = FALSE)

write.csv(pretty_table,
          "data/power.csv",
          row.names = FALSE)

##Thanks nice table
#MSE
#$MCSE_MSE
table_bias_empse <- perf_clean %>%
  filter(n %in% ns_to_show) %>% # ---- FIX METHOD LABEL HERE ----
mutate(
  Method = case_when(
    Method == "CITS"              ~ "CITS",
    Method == "CITS_0.4_constant" ~ "CITS(rho = 0.4)",
    Method == "Trd"               ~ "Trd",
    TRUE                          ~ as.character(Method)
  )
)%>%
  
  # string for each cell: Estimate (Empirical_SE)
  mutate(
    cell = sprintf("%.4f (%.4f)", MSE_estimate, MCSE_MSE),
    
    # column key = Effect-with-true-val + rho
    col_key = paste0(effect_label, "_rho", rho)
  ) %>%
  
  select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small(-0.0408)", "Moderate(-0.3567)", "Large(-0.5108)")

col_order <- c(
  "n", "Method",
  unlist(lapply(effects, function(eff) paste0(eff, "_rho", rho_seq)))
)

estimates <- table_bias_empse[, col_order]

# Row 1: effect names (repeated for each rho)
`Effect size`<- c(
  "Length of series (n)",          # col 1
  "Method",                        # col 2
  rep("Small effect (-0.0408)",    length(rho_seq)),
  rep("Moderate effect (-0.3567)", length(rho_seq)),
  rep("Large effect (-0.5108)",    length(rho_seq))
)

# Row 2: rho values under each effect
`Correlation`<- c(
  "",                              # under n
  "",                              # under Method
  rep(rho_seq, times = length(effects))
)

# Coerce everything to character so we can rbind
estimates_chr <- estimates %>% mutate(across(everything(), as.character))

pretty_table <- rbind(
  `Effect size`,
  `Correlation`,
  as.matrix(estimates_chr)
)

# Optional: remove column names so Excel doesn’t show them
colnames(pretty_table) <- rep("", ncol(pretty_table))

View(pretty_table)
write.csv(pretty_table,
          "data/mse.csv",
          row.names = FALSE)

write.csv(pretty_table,
          "data/mse.csv",
          row.names = FALSE)

##Emperical and b=model based SE
table_bias_empse <- perf_clean %>%
  filter(n %in% ns_to_show) %>% # ---- FIX METHOD LABEL HERE ----
mutate(
  Method = case_when(
    Method == "CITS"              ~ "CITS",
    Method == "CITS_0.4_constant" ~ "CITS(rho = 0.4)",
    Method == "Trd"               ~ "Trd",
    TRUE                          ~ as.character(Method)
  )
)%>%
  
  # string for each cell: Estimate (Empirical_SE)
  mutate(
    cell = sprintf("%.4f (%.4f)", Empirical_SE, avg_model_SE),
    
    # column key = Effect-with-true-val + rho
    col_key = paste0(effect_label, "_rho", rho)
  ) %>%
  
  select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small(-0.0408)", "Moderate(-0.3567)", "Large(-0.5108)")

col_order <- c(
  "n", "Method",
  unlist(lapply(effects, function(eff) paste0(eff, "_rho", rho_seq)))
)

estimates <- table_bias_empse[, col_order]

# Row 1: effect names (repeated for each rho)
`Effect size`<- c(
  "Length of series (n)",          # col 1
  "Method",                        # col 2
  rep("Small effect (-0.0408)",    length(rho_seq)),
  rep("Moderate effect (-0.3567)", length(rho_seq)),
  rep("Large effect (-0.5108)",    length(rho_seq))
)

# Row 2: rho values under each effect
`Correlation`<- c(
  "",                              # under n
  "",                              # under Method
  rep(rho_seq, times = length(effects))
)

# Coerce everything to character so we can rbind
estimates_chr <- estimates %>% mutate(across(everything(), as.character))

pretty_table <- rbind(
  `Effect size`,
  `Correlation`,
  as.matrix(estimates_chr)
)

# Optional: remove column names so Excel doesn’t show them
colnames(pretty_table) <- rep("", ncol(pretty_table))

write.csv(pretty_table,
          "data/SE.csv",
          row.names = FALSE)

write.csv(pretty_table,
          "data/SE.csv",
          row.names = FALSE)



