
#packages
library(dplyr)
library(cowplot)
library(patchwork)
#library(ggplot2)
library(grid)
library(ggplot2)
library(tidyr)
library(stringr)
if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}
library(ggridges)
##OBJECTIVE 1
#perf<-readRDS("data/performance_summary.rds")
perf <- read.csv("data/performance_summary.csv")
str(perf)


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
      factor(Method, levels = c("Trd", "CITS","CITS_0.4_constant","CITS_spill_3pct"))
    } else {
      Method
    },
    
    ## --- dynamic facet labels using true_val ---
    facet_label = case_when(
      EstimandScenario == "Small"    ~ paste0("Small (", round(true_val, 4), ")"),
      EstimandScenario == "Moderate" ~ paste0("Moderate (", round(true_val, 4), ")"),
      EstimandScenario == "Large"    ~ paste0("Large (", round(true_val, 4), ")"),
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

#Empirical and model based standard errors
##traditional regression plot
p_trd2 <- perf_plot %>%
  filter(Method == "CITS_spill_3pct") %>%
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

p3<-p_trd2

ggsave(
  filename = "graphs/p_cits_SE_10pctspill.png",
  plot     = p_trd2,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)

ggsave(
  filename = "graphs/p_cits_SE_10pctspill.pdf",
  plot     = p_trd2,
  width    = 8,      # inches
  height   = 4.5,
  device   = cairo_pdf
)




p_trd2 <- perf_plot %>%
  filter(Method == "CITS_spill_3pct") %>%
  ggplot(aes(x = n, group = 1)) +
  
  geom_line(
    aes(y = Empirical_SE, colour = "Empirical SE"),
    linewidth = 0.7
  ) +
  
  geom_point(
    aes(y = Empirical_SE, colour = "Empirical SE"),
    size = 1.8
  ) +
  
  geom_line(
    aes(y = avg_model_SE, colour = "Model-based SE"),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_point(
    aes(y = avg_model_SE, colour = "Model-based SE"),
    size = 1.8
  ) +
  
  facet_grid(facet_label ~ rho_f) +
  
  scale_colour_manual(
    name = "SE type",
    values = c(
      "Empirical SE"   = "#F8766D",
      "Model-based SE" = "#00BFC4"
    )
  ) +
  
  labs(
    x = "Length of time series (n)",
    y = "Standard error (SE)",
    title = ""
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    # Top strips: autocorrelation levels
    strip.text.x = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 4, b = 4)
    ),
    
    # Right strips: Small, Moderate, Large
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t =3, b =3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p_trd2

# Save with slightly larger width (important)
ggsave(
  filename = "graphs/p_cits_SE_10pctspill_1.png",
  plot     = p_trd2,
  width    = 11,
  height   = 5.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_cits_SE_10pctspill_1.pdf",
  plot     = p_trd2,
  width    = 11,
  height   = 5.5,
  device   = cairo_pdf
)





##Coverage graphs
##Confidence interval coverage
#coveragetrd
p_trd_cova <- perf_plot %>%
  filter(Method == "CITS_spill_3pct") %>%
  ggplot(
    aes(
      x      = n,
      y      = coverage_pct,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  geom_hline(
    yintercept = 95,
    colour = "black",
    linetype = "solid",
    linewidth = 0.8
  ) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Coverage (%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p5 <- p_trd_cova

ggsave(
  filename = "graphs/p_CITS_spill_coverage.png",
  plot     =p_trd_cova ,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)


ggsave(
  filename = "graphs/p_CITS_spill_coverage.pdf",
  plot     =p_trd_cova ,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  device =cairo_pdf     # good for publication
)




#library(dplyr)

coverage_summary <- perf %>%
  filter(Method == "CITS_spill_3pct") %>%
  group_by(EstimandScenario) %>%
  summarise(
    mean_coverage = mean(coverage_pct, na.rm = TRUE),
    sd_coverage   = sd(coverage_pct, na.rm = TRUE),
    min_coverage  = min(coverage_pct, na.rm = TRUE),
    max_coverage  = max(coverage_pct, na.rm = TRUE),
    .groups = "drop"
  )

coverage_summary



###Power plot Multivariable regression
##traditional regression plot
#rho_f<-c(0,0.2,0.4,0.6,0.8)
p_trd_power <- perf_plot %>%
  filter(Method == "CITS_spill_3pct") %>%
  ggplot(
    aes(
      x      = n,
      y      = power_pct,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  geom_hline(
    yintercept = 80,
    linetype = "solid",
    colour = "black",
    linewidth = 0.8
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Power (%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p7 <- p_trd_power
p7
ggsave(
  filename = "graphs/p_CITS_spill_Power.png",
  plot     =p7,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)


ggsave(
  filename = "graphs/p_CITS_spill_Power.pdf",
  plot     =p7,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  device=cairo_pdf    # good for publication
)






#MSE
##trd MSE
p_trd_mse <- perf_plot %>%
  filter(Method == "CITS_spill_3pct") %>%
  ggplot(
    aes(
      x      = n,
      y      = MSE_estimate,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  # CHANGED: slightly thicker lines for consistency with the power plot
  # and better readability in the thesis PDF.
  geom_line(
    linewidth = 0.7
  ) +
  
  # CHANGED: slightly larger points so the individual simulation
  # conditions remain visible when the figure is reduced on the page.
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    
    # CHANGED: corrected spacing and capitalization.
    y = "Mean squared error (MSE)",
    
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  
  # ADDED: use the same rho notation as in the power plot
  # so all simulation figures are visually consistent.
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  # ADDED: matching shape labels for consistency.
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    # CHANGED: match the facet-strip formatting used for your power plot.
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    # ADDED: explicit font sizes so the figure remains readable
    # after insertion into the thesis.
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p9 <- p_trd_mse
p9

ggsave(
  filename = "graphs/p_CITS_spill_mse.png",
  plot     =p_trd_mse,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)


ggsave(
  filename = "graphs/p_CITS_spill_mse.pdf",
  plot     =p_trd_mse,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  device      =cairo_pdf     # good for publication
)



##Plotting point estimates
## If you saved estimates earlier:
## 1. Base data (what you already have) write.csv(perf, file = "data/perf.csv", row.names = FALSE)
#estimates <- readRDS("data/sim_estimates.rds")
estimates<- read.csv("data/sim_estimates.csv")
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
      levels = c("0", "0.2", "0.4", "0.6", "0.8")
      #labels = c("ρ=0.0", "ρ=0.2", "ρ=0.4", "ρ=0.6", "ρ=0.8")
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
  filter(Method == "CITS_spill_3pct", n %in% selected_n) %>%
  left_join(
    scenario_labels %>%dplyr::select(EstimandScenario, facet_label),
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
    
    # CHANGED: use a thinner outline because five rho distributions
    # overlap within each time-series length. This reduces visual crowding.
    linewidth      = 0.6,
    
    # CHANGED: slightly more transparency so overlapping ridgelines
    # are easier to distinguish.
    alpha          = 0.6
  ) +
  
  geom_vline(
    data        = true_lines,
    aes(xintercept = true_B),
    linetype    = "dashed",
    colour      = "black",
    
    # ADDED: slightly stronger line width so the true effect remains
    # clearly visible when printed.
    linewidth   = 0.7
    
    #inherit.aes = FALSE
  ) +
  
  facet_grid(. ~ facet_label) +
  
  labs(
    title = "",
    
    x = "Estimated intervention effect (log IRR)",
    
    # CHANGED: shorter and slightly more direct wording.
    # n still makes clear that this is the series length.
    y = "Time-series length (n)"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  ggridges::theme_ridges() +
  
  theme(
    legend.position = "bottom",
    
    strip.background = element_rect(fill = "grey95"),
    
    # CHANGED: slightly larger facet labels for better readability
    # in the printed thesis.
    strip.text = element_text(
      size = 10,
      face = "bold",
      margin = margin(t = 3, b = 4)
    ),
    
    # ADDED: explicitly increase axis-title size rather than relying
    # on theme_ridges() defaults.
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    # ADDED: enlarge tick labels so values remain readable when the
    # figure is inserted at thesis-page width.
    axis.text = element_text(
      size = 11
    ),
    
    # ADDED: improve legend readability in print.
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p_trd_ridge
p1<-p_trd_ridge
ggsave(
  "graphs/CITS_spill_ridgeline_estimates.png",
  p_trd_ridge,
  width  = 8,
  height = 4.5,
  dpi    = 300
)

ggsave(
  filename = "graphs/CITS_spill_ridgeline_estimates.pdf",
  plot     = p_trd_ridge,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)








##Bias
p_trd_bias <- perf_plot %>%
  filter(Method == "CITS_spill_3pct") %>%   # CHECK: should this be 10pct?
  ggplot(
    aes(
      x      = n,
      y      = Bias_percent,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    title = "",
    x = "Length of time series (n)",
    y = "Percent bias (%)"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 5, b = 5, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p12 <- p_trd_bias

ggsave(
  filename = "graphs/p_CITS_spill_bias.png",
  plot     = p_trd_bias,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)

ggsave(
  filename = "graphs/p_CITS_spill_bias.pdf",
  plot     = p_trd_bias,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  device   = cairo_pdf    # good for publication
)




Bias_summary <- perf %>%
  filter(Method == "CITS_spill_3pct") %>%
  group_by(EstimandScenario,rho) %>%
  summarise(
    mean_bias= mean(Bias_percent, na.rm = TRUE),
    sd_bias= sd(Bias_percent, na.rm = TRUE),
    min_bias= min(Bias_percent, na.rm = TRUE),
    max_bias= max(Bias_percent, na.rm = TRUE),
    .groups = "drop"
  )

Bias_summary





##lets create some table 
ns_to_show<-c(c(12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,80,100,150))
perf<-perf%>%filter(Method!="CITS_spill_3pct")
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
      EstimandScenario == "Small"    ~ paste0("Small (",    round(true_val, 4), ")"),
      EstimandScenario == "Moderate" ~ paste0("Moderate (", round(true_val, 4), ")"),
      EstimandScenario == "Large"    ~ paste0("Large (",    round(true_val, 4), ")"),
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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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









##OBJECTIVE 2
##reading data
#perf<-readRDS("data/performance_summary.rds")
perf <- read.csv("data/performance_summary.csv")
str(perf)

#perf<-filter(perf,Method!="CITS_spill_3pct")
##creating a table of performance measures to be used in R_shinydashbord
perfRshiny <- perf %>%
  mutate(
    `Policy effect (Empirical SE)`        = sprintf("%.4f (%.4f)", B_hat, Empirical_SE),
    `Bias (MCSE)`                         = sprintf("%.4f (%.4f)", Bias, MCSE_Bias),
    `Coverage 95% (MCSE)`                 = sprintf("%.2f (%.4f)", coverage_pct, MCSE_cov_pct),
    `Power (MCSE)`                        = sprintf("%.2f (%.4f)", power_pct, MCSE_power_pct),
    `Empirical SE (Model-based SE)`       = sprintf("%.4f (%.4f)", Empirical_SE, avg_model_SE),
    `Mean square error (MCSE)`            = sprintf("%.4f (%.4f)", MSE_estimate, MCSE_MSE),
    `Percent Bias`=Bias_percent,
    Method = case_when(
      Method == "CITS_0.4_constant" ~ "Controlled interrupted time series(CITS,ρ = 0.4)",
      Method == "Trd"               ~ "Multivariable negative binomial regression(Trd)",
      Method == "CITS_spill_3pct"               ~ "Controlled interrupted time series (Contaminated control)",
      TRUE                          ~ "Controlled interrupted time series(CITS)"
    )
  ) %>%dplyr::select(
    rho, n, EstimandScenario, Method,
    `Policy effect (Empirical SE)`,
    `Bias (MCSE)`,
    `Coverage 95% (MCSE)`,
    `Power (MCSE)`,
    `Empirical SE (Model-based SE)`,
    `Mean square error (MCSE)`,
    true_val,
    `Percent Bias`
  )

table(perfRshiny$Method,useNA ="ifany")
perfRshiny<-perfRshiny%>%filter(Method!="Controlled interrupted time series (Contaminated control)")



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
      factor(Method, levels = c("Trd", "CITS","CITS_0.4_constant","CITS_spill_3pct"))
    } else {
      Method
    },
    
    ## --- dynamic facet labels using true_val ---
    facet_label = case_when(
      EstimandScenario == "Small"    ~ paste0("Small (", round(true_val, 4), ")"),
      EstimandScenario == "Moderate" ~ paste0("Moderate (", round(true_val, 4), ")"),
      EstimandScenario == "Large"    ~ paste0("Large (", round(true_val, 4), ")"),
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
geom_line(
  aes(y = Empirical_SE, colour = "Empirical SE"),
  linewidth = 0.7
) +
  
  geom_point(
    aes(y = Empirical_SE, colour = "Empirical SE"),
    size = 1.8
  ) +
  
  # ---- Model-based SE ----
geom_line(
  aes(y = avg_model_SE, colour = "Model-based SE"),
  linetype = "dashed",
  linewidth = 0.7
) +
  
  geom_point(
    aes(y = avg_model_SE, colour = "Model-based SE"),
    size = 1.8
  ) +
  
  facet_grid(facet_label ~ rho_f) +
  
  scale_colour_manual(
    name = "SE type",
    values = c(
      "Empirical SE"   = "#F8766D",
      "Model-based SE" = "#00BFC4"
    )
  ) +
  
  labs(
    x = "Length of time series (n)",
    y = "Standard error (SE)",
    title = "A"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.x = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 4, b = 4)
    ),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p3 <- p_trd2

ggsave(
  filename = "graphs/p_trd_SE_traditional.png",
  plot     = p_trd2,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_trd_SE_traditional.pdf",
  plot     = p_trd2,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


##CITS Controlled Interrupted Time Series (CITS)
p_cits2 <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(aes(x = n, group = 1)) +
  
  # ---- Empirical SE ----
geom_line(
  aes(y = Empirical_SE, colour = "Empirical SE"),
  linewidth = 0.7
) +
  
  geom_point(
    aes(y = Empirical_SE, colour = "Empirical SE"),
    size = 1.8
  ) +
  
  # ---- Model-based SE ----
geom_line(
  aes(y = avg_model_SE, colour = "Model-based SE"),
  linetype = "dashed",
  linewidth = 0.7
) +
  
  geom_point(
    aes(y = avg_model_SE, colour = "Model-based SE"),
    size = 1.8
  ) +
  
  facet_grid(facet_label ~ rho_f) +
  
  scale_colour_manual(
    name = "SE type",
    values = c(
      "Empirical SE"   = "#F8766D",
      "Model-based SE" = "#00BFC4"
    )
  ) +
  
  labs(
    x = "Length of time series (n)",
    y = "Standard error (SE)",
    title = "B"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.x = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 4, b = 4)
    ),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.position = "none",
    
    panel.grid.minor = element_blank()
  )

p4 <- p_cits2

ggsave(
  filename = "graphs/p_cits_SE.png",
  plot     = p_cits2,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_cits_SE.pdf",
  plot     = p_cits2,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


combinedstd <- p3 / p4

ggsave(
  filename = "graphs/SE_both.pdf",
  plot     = combinedstd,
  width    = 11,
  height   = 12,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "graphs/SE_both.png",
  plot     = combinedstd,
  width    = 11,
  height   = 12,
  units    = "in",
  dpi      = 300
)


 
#CITS with fixed autocorelation(0.4)
##CITS Controlled Interrupted Time Series (CITS)
table(perf_plot$Method)
p_cits2 <- perf_plot %>%
  filter(Method == "CITS_0.4_constant") %>%
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
    title = "A"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )

p_cits2
ggsave(
  filename = "graphs/p_cits0.4_SE.png",
  plot     = p_cits2,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)

ggsave(
  filename = "graphs/p_cits0.4_SE.pdf",
  plot     = p_cits2,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)







#original scale
#log scale
table(perf_plot$Method, useNA = "ifany")

# Keep only multivariable regression and CITS
perf_plotcitstrd <- perf_plot %>%
  filter(Method == "Trd" | Method == "CITS")

p_combined_log2 <- ggplot(
  perf_plotcitstrd,
  aes(
    x = n,
    y = Empirical_SE,
    colour = Method,
    group = Method
  )
) +
  
  # CHANGED: slightly thicker lines for better visibility
  # in the thesis and during presentation
  geom_line(
    linewidth = 0.7
  ) +
  
  # CHANGED: slightly larger points for readability
  geom_point(
    size = 2
  ) +
  
  facet_grid(facet_label ~ rho_f) +
  
  # Retain log scale because empirical SE varies substantially
  # across simulation conditions
  scale_y_log10() +
  
  scale_colour_discrete(
    labels = c(
      Trd = "Multivariable NB",
      CITS = "CITS"
    )
  ) +
  
  labs(
    x = "Length of time series (n)",
    
    # CHANGED: improved spacing/formatting
    y = "Empirical SE (log scale)",
    
    colour = "Method"
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    # ADDED: improve readability of facet labels
    strip.text = element_text(
      size = 9,
      face = "bold"
    ),
    
    # ADDED: improve axis-title readability
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    # ADDED: larger axis values
    axis.text = element_text(
      size = 10
    ),
    
    # ADDED: improve legend readability
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.position = "right",
    
    panel.grid.minor = element_blank()
  )

p_combined_log2


ggsave(
  filename = "graphs/p_trd_SE_cits_trdCombined.png",
  plot     = p_combined_log2,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)


ggsave(
  filename = "graphs/p_trd_SE_cits_trdCombined.pdf",
  plot     = p_combined_log2,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)
##All cits
#perf_plotcitstrd<-perf_plot[perf_plot$Method!="CITS_0.4_constant",]
perf_plot11<-perf_plot%>%filter(Method!="CITS_spill_3pct")
perf_plot1<-perf_plot11%>%mutate(Method = dplyr::case_when(
  Method == "CITS_0.4_constant" ~ "CITS (ρ=0.4)",
  Method == "CITS"             ~ "CITS",
  Method == "Trd"              ~ "Trd",
  TRUE                         ~ Method
))
p_combined_log2 <- ggplot(
  perf_plot1,
  aes(x = n, y = Empirical_SE, colour = Method, group = Method)
) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(facet_label ~ rho_f) +
  scale_y_log10() +
  scale_colour_discrete(labels = c(Trd = "Multivariable NB", CITS = "CITS")) +
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
  filename = "graphs/p_trd_SE_citstypes_trdCombined.png",
  plot     = p_combined_log2,
  width    = 8,      # inches
  height   = 4.5,    # adjust as you like
  dpi      = 300     # good for publication
)


ggsave(
  filename = "graphs/p_trd_SE_citstypes_trdCombined.pdf",
  plot     = p_combined_log2,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)





##Coverage graphs
#coveragetrd
p_trd_cova <- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(
    aes(
      x      = n,
      y      = coverage_pct,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  # 95% nominal coverage reference line
  geom_hline(
    yintercept = 95,
    colour = "black",
    linetype = "solid",
    linewidth = 0.8
  ) +
  
  # Slightly thicker lines for readability
  geom_line(
    linewidth = 0.7
  ) +
  
  # Slightly larger points
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Coverage (%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "A"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p5 <- p_trd_cova

ggsave(
  filename = "graphs/p_trd_coverage.png",
  plot     = p_trd_cova,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_trd_coverage.pdf",
  plot     = p_trd_cova,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


# CITS
p_trd_cits <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      = coverage_pct,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  geom_hline(
    yintercept = 95,
    colour = "black",
    linetype = "solid",
    linewidth = 0.8
  ) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Coverage (%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "B"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "none",
    
    panel.grid.minor = element_blank()
  )

p6 <- p_trd_cits

ggsave(
  filename = "graphs/p_cits_coverage.png",
  plot     = p_trd_cits,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_cits_coverage.pdf",
  plot     = p_trd_cits,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


combinedcova <- p5 / p6

ggsave(
  filename = "graphs/cova_both.pdf",
  plot     = combinedcova,
  width    = 11,
  height   = 12,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "graphs/cova_both.png",
  plot     = combinedcova,
  width    = 11,
  height   = 12,
  units    = "in",
  dpi      = 300
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
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  geom_hline(
    yintercept = 80,
    linetype = "solid",
    colour = "black",
    linewidth = 0.8
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Power (%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "A"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p7 <- p_trd_power

ggsave(
  filename = "graphs/p_trd_Power.png",
  plot     = p_trd_power,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_trd_Power.pdf",
  plot     = p_trd_power,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


# CITS
p_cits_power <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      = power_pct,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  geom_hline(
    yintercept = 80,
    linetype = "solid",
    colour = "black",
    linewidth = 0.8
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Power (%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "B"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "none",
    
    panel.grid.minor = element_blank()
  )

p8 <- p_cits_power

ggsave(
  filename = "graphs/p_cits_power.png",
  plot     = p_cits_power,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_cits_power.pdf",
  plot     = p_cits_power,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


combinedpowa <- p7 / p8

ggsave(
  filename = "graphs/powa_both.pdf",
  plot     = combinedpowa,
  width    = 11,
  height   = 12,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "graphs/powa_both.png",
  plot     = combinedpowa,
  width    = 11,
  height   = 12,
  units    = "in",
  dpi      = 300
)


#MSE
##trd MSE
p_trd_mse <- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(
    aes(
      x      = n,
      y      = MSE_estimate,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  # Slightly stronger lines for readability and consistency
  geom_line(
    linewidth = 0.7
  ) +
  
  # Slightly larger points for thesis/presentation readability
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Mean squared error (MSE)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "A"
  ) +
  
  # Consistent rho notation across all Objective 2 figures
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p9 <- p_trd_mse


ggsave(
  filename = "graphs/p_trd_mse.png",
  plot     = p_trd_mse,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)


ggsave(
  filename = "graphs/p_trd_mse.pdf",
  plot     = p_trd_mse,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


# CITS
p_cits_mse <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      = MSE_estimate,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  # Same line width as Panel A
  geom_line(
    linewidth = 0.7
  ) +
  
  # Same point size as Panel A
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    y = "Mean squared error (MSE)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "B"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 8,
      face = "bold",
      margin = margin(t = 3, b = 3, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    # Legend shown only in Panel A
    legend.position = "none",
    
    panel.grid.minor = element_blank()
  )

p10 <- p_cits_mse


ggsave(
  filename = "graphs/p_cits_mse.png",
  plot     = p_cits_mse,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)


ggsave(
  filename = "graphs/p_cits_mse.pdf",
  plot     = p_cits_mse,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


combinedmse <- p9 / p10


ggsave(
  filename = "graphs/mse_both.pdf",
  plot     = combinedmse,
  width    = 12,
  height   = 14,
  units    = "in",
  device   = cairo_pdf
)


ggsave(
  filename = "graphs/mse_both.png",
  plot     = combinedmse,
  width    = 12,
  height   = 14,
  units    = "in",
  dpi      = 300
)

##Plotting point estimates
## If you saved estimates earlier:
## 1. Base data (what you already have) write.csv(perf, file = "data/perf.csv", row.names = FALSE)
#estimates <- readRDS("data/sim_estimates.rds")
estimates<- read.csv("data/sim_estimates.csv")
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
      levels = c("0", "0.2", "0.4", "0.6", "0.8")
      #labels = c("ρ=0.0", "ρ=0.2", "ρ=0.4", "ρ=0.6", "ρ=0.8")
    ),
    n = if (!is.numeric(n)) as.numeric(as.character(n)) else n
  ) %>%
  dplyr::select(-rho_char)

## 2. Scenario labels + true effects (fix facet order here)
scenario_labels <- est_plot %>%
  distinct(EstimandScenario, true_B) %>%
  arrange(EstimandScenario) %>%
  mutate(
    facet_label = paste0(
      EstimandScenario, " (", round(true_B, 4), ")"
    ),
    facet_label = factor(facet_label, levels = facet_label)
  )


## This will also be used for the dashed lines
true_lines <- scenario_labels

## 3. Subset to selected n and join facet labels
selected_n <- c(12, 24, 40, 60, 80, 100)

est_trd <- est_plot %>%
  filter(Method == "Trd", n %in% selected_n) %>%
  left_join(
    scenario_labels %>%dplyr::select(EstimandScenario, facet_label),
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
    
    # CHANGED: thinner outline to reduce crowding
    linewidth      = 0.6,
    
    # CHANGED: slightly more transparent for overlapping densities
    alpha          = 0.6
  ) +
  
  geom_vline(
    data        = true_lines,
    aes(xintercept = true_B),
    linetype    = "dashed",
    colour      = "black",
    
    # CHANGED: clearer true-effect reference line
    linewidth   = 0.7
  ) +
  
  facet_grid(. ~ facet_label) +
  
  # ADDED: same x-axis range as CITS for fair visual comparison
  coord_cartesian(xlim = c(-2, 2)) +
  
  labs(
    title = "A",
    x     = NULL,
    y     = "Time-series length (n)"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  ggridges::theme_ridges() +
  
  theme(
    legend.position = "none",
    
    strip.background = element_rect(fill = "grey95"),
    
    # CHANGED: slightly cleaner facet labels
    strip.text = element_text(
      size = 10,
      face = "bold",
      margin = margin(t = 3, b = 4)
    ),
    
    # ADDED: prevents letters such as g from being clipped
    strip.clip = "off",
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    panel.grid.minor = element_blank()
  )

p_trd_ridge
p1 <- p_trd_ridge

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
      levels = c("0", "0.2", "0.4", "0.6", "0.8")
      #labels = c("ρ = 0.0", "ρ = 0.2", "ρ = 0.4", "ρ = 0.6", "ρ = 0.8")
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
  filter(Method == "CITS", n %in% selected_n) %>%
  left_join(
    scenario_labels %>% dplyr::select(EstimandScenario, facet_label),
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
    
    # CHANGED: thinner outline for consistency with Panel A
    linewidth      = 0.6,
    
    # CHANGED: same transparency as Panel A
    alpha          = 0.6
  ) +
  
  geom_vline(
    data        = true_lines,
    aes(xintercept = true_B),
    linetype    = "dashed",
    colour      = "black",
    
    # CHANGED: clearer true-effect line
    linewidth   = 0.7
  ) +
  
  facet_grid(. ~ facet_label) +
  
  # ADDED: same x-axis range as Panel A
  coord_cartesian(xlim = c(-2, 2)) +
  
  labs(
    title = "B",
    x     = "Estimated intervention effect (log IRR)",
    y     = "Time-series length (n)"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  ggridges::theme_ridges() +
  
  theme(
    legend.position = "bottom",
    
    strip.background = element_rect(fill = "grey95"),
    
    strip.text = element_text(
      size = 10,
      face = "bold",
      margin = margin(t = 3, b = 4)
    ),
    
    strip.clip = "off",
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p_trd_ridge
p2 <- p_trd_ridge

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
  device = cairo_pdf
)

##lets try doing graph combining the ridgeline_plots
# p1 and p2 are your already-finished plots
combined <- p1 / p2

ggsave(
  "graphs/ridgeline_both.pdf",
  plot = combined,
  width = 12,
  height = 14,
  units = "in",
  dpi = 300
)

ggsave(
  "graphs/ridgeline_both.png",
  plot = combined,
  width = 12,
  height = 14,
  units = "in",
  dpi = 300
)

combined2 <- plot_grid(
  p1, p2,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 14,
  align = "v"
)

ggsave(
  "ridgeline_bothCow.pdf",
  plot = combined2,
  width = 12,
  height = 14,
  units = "in",
  dpi = 300
)




#CITS with autocorelation fixe at 0.4
table(est_plot$Method)
est_plot1<-est_plot%>%mutate(Method = case_when(
  Method == "CITS"              ~ "CITS",
  Method == "CITS_0.4_constant" ~ "CITS(ρ= 0.4)",
  Method == "Trd"               ~ "Trd",
  TRUE                          ~ as.character(Method)
))
table(est_plot1$Method)
est_trd <- est_plot1 %>%
  filter(Method == "CITS(ρ= 0.4)", n %in% selected_n) %>%
  left_join(
    scenario_labels %>%dplyr::select(EstimandScenario, facet_label),
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
  "graphs/CITS0.4_ridgeline_estimates.png",
  p_trd_ridge,
  width  = 8,
  height = 4.5,
  dpi    = 300
)

ggsave(
  "graphs/CITS0.4_ridgeline_estimates.pdf",
  p_trd_ridge,
  width  = 8,
  height = 4.5,
  device   = cairo_pdf
)







##Bias
p_trd_bias <- perf_plot %>%
  filter(Method == "Trd") %>%
  ggplot(
    aes(
      x      = n,
      y      = Bias_percent,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +# CHANGED: slightly stronger lines for thesis readability
  geom_line(
    linewidth = 0.7
  ) +
  
  # CHANGED: slightly larger points
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  # ADDED: use the same y-axis range for the methods being compared
  # so visual comparisons between A and B are fair
  coord_cartesian(
    ylim = range(
      perf_plot$Bias_percent[
        perf_plot$Method %in% c("Trd", "CITS")
      ],
      na.rm = TRUE
    )
  ) +
  
  labs(
    x = "Length of time series (n)",
    
    # CHANGED: clearer statistical terminology
    y = "Percent bias (%)",
    
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "A"
  ) +
  
  # ADDED: consistent rho notation across all your simulation plots
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 5, b = 5, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    panel.grid.minor = element_blank()
  )

p12 <- p_trd_bias

ggsave(
  filename = "graphs/p_trd_bias.png",
  plot     = p_trd_bias,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_trd_bias.pdf",
  plot     = p_trd_bias,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)



#cits
p_cits_bias <- perf_plot %>%
  filter(Method == "CITS") %>%
  ggplot(
    aes(
      x      = n,
      y      = Bias_percent,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) + geom_line(
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  # SAME y-axis range as Panel A
  coord_cartesian(
    ylim = range(
      perf_plot$Bias_percent[
        perf_plot$Method %in% c("Trd", "CITS")
      ],
      na.rm = TRUE
    )
  ) +
  
  labs(
    x = "Length of time series (n)",
    y = "Percent bias (%)",
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "B"
  ) +
  
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    strip.text.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 5, b = 5, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    # Keep no legend here because Panel A already contains it
    legend.position = "none",
    
    panel.grid.minor = element_blank()
  )

p13 <- p_cits_bias

ggsave(
  filename = "graphs/p_cits_bias.png",
  plot     = p_cits_bias,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_cits_bias.pdf",
  plot     = p_cits_bias,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)

##CITS with constant corelation
table(perf_plot$Method)

p_cits_bias <- perf_plot %>%
  filter(Method == "CITS_0.4_constant") %>%
  ggplot(
    aes(
      x      = n,
      y      = Bias_percent,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  # CHANGED: slightly thicker lines for consistency
  geom_line(
    linewidth = 0.7
  ) +
  
  # CHANGED: slightly larger points for readability
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    
    # CHANGED: clearer wording
    y = "Percent bias (%)",
    
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = "C"
  ) +
  
  # ADDED: consistent rho labels
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  # ADDED: matching shape labels
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    # CHANGED: improve facet-label readability
    strip.text.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 5, b = 5, r = 6, l = 6)
    ),
    
    # ADDED: explicit axis-title size
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    # ADDED: larger tick labels
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "none",
    
    panel.grid.minor = element_blank()
  )

p14 <- p_cits_bias

ggsave(
  filename = "graphs/p_citsrho0.4_bias.png",
  plot     = p_cits_bias,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_citsrho0.4_bias.pdf",
  plot     = p_cits_bias,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)


p_cits_bias1 <- perf_plot %>%
  filter(Method == "CITS_0.4_constant") %>%
  ggplot(
    aes(
      x      = n,
      y      = Bias_percent,
      colour = rho_f,
      shape  = rho_f,
      group  = rho_f
    )
  ) +
  
  # CHANGED: same line thickness as the other bias plots
  geom_line(
    linewidth = 0.7
  ) +
  
  # CHANGED: same point size as the other bias plots
  geom_point(
    size = 2.2
  ) +
  
  facet_grid(facet_label ~ .) +
  
  labs(
    x = "Length of time series (n)",
    
    # CHANGED: clearer wording
    y = "Percent bias (%)",
    
    colour = "Autocorrelation (ρ)",
    shape  = "Autocorrelation (ρ)",
    title = ""
  ) +
  
  # ADDED: consistent rho labels
  scale_colour_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  # ADDED: matching shape labels
  scale_shape_discrete(
    labels = c(
      expression(rho == 0.0),
      expression(rho == 0.2),
      expression(rho == 0.4),
      expression(rho == 0.6),
      expression(rho == 0.8)
    ),
    name = expression("Autocorrelation (" * rho * ")")
  ) +
  
  theme_bw() +
  
  theme(
    strip.background = element_rect(fill = "grey95"),
    
    # CHANGED: same facet-strip formatting as above
    strip.text.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 5, b = 5, r = 6, l = 6)
    ),
    
    axis.title.x = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    legend.position = "none",
    
    panel.grid.minor = element_blank()
  )


ggsave(
  filename = "graphs/p_citsrho0.4_bias1.png",
  plot     = p_cits_bias1,
  width    = 8,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = "graphs/p_citsrho0.4_bias1.pdf",
  plot     = p_cits_bias1,
  width    = 8,
  height   = 4.5,
  device   = cairo_pdf
)

#combine
combinedBiased<- p12 / p13 / p14
ggsave(
  "graphs/combinedBIAS_plots.pdf",
  plot   = combinedBiased,
  width  = 12,
  height = 14,
  units  = "in",
  device = cairo_pdf
)

ggsave(
  "graphs/combinedBIAS_plots.png",
  plot   = combinedBiased,
  width  = 12,
  height = 14,
  units  = "in",
  dpi    = 300
)

combinedBiasedThesis<- p12 / p13


ggsave(
  "graphs/combinedBIASthesis_plots.pdf",
  plot   =combinedBiasedThesis,
  width  = 11,
  height =12,
  units  = "in",
  device = cairo_pdf
)

ggsave(
  "graphs/combinedBIASthesis_plots.png",
  plot   =combinedBiasedThesis,
  width  = 11,
  height =12,
  units  = "in",
  dpi    = 300
)



##lets create some table 
ns_to_show<-c(c(12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,80,100,150))
perf<-perf%>%filter(Method!="CITS_spill_3pct")
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
      EstimandScenario == "Small"    ~ paste0("Small (",    round(true_val, 4), ")"),
      EstimandScenario == "Moderate" ~ paste0("Moderate (", round(true_val, 4), ")"),
      EstimandScenario == "Large"    ~ paste0("Large (",    round(true_val, 4), ")"),
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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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
  
  dplyr::select(n, Method, col_key, cell) %>%
  distinct() %>%                      # just in case
  
  pivot_wider(
    names_from  = col_key,
    values_from = cell
  ) %>%
  arrange(n, Method)

table_bias_empse


rho_seq    <- c("0.0", "0.2", "0.4", "0.6", "0.8")
effects    <- c("Small (-0.0408)", "Moderate (-0.3567)", "Large (-0.5108)")

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



##OBjective1 thesis
 ##Bias of estimated intervebntion effect(percent bias and ridgeline plot)
est_trd <- est_plot %>%filter(Method == "CITS_spill_3pct")
est_trd<-est_trd%>%arrange(true_B,rho,n)%>%select(true_B,estimate,rho,n,everything())
#library(dplyr)

bias_dat<-est_trd

bias_dat <-bias_dat%>%
  mutate(
    bias = estimate - true_B,
    abs_bias = abs(bias),
    pct_bias = 100 * abs(bias / true_B)
  )

abs_bias_effect <- bias_dat %>%
  group_by(EstimandScenario) %>%
  summarise(
    mean_abs_bias =round(mean(abs_bias, na.rm = TRUE),2),
    sd_abs_bias =round(sd(abs_bias, na.rm = TRUE),2),
    mean_estimate=round(mean(estimate,na.rm=TRUE),2),
    .groups = "drop"
  )

abs_bias_effect


p_trd_biasNHJ<- perf_plot %>%filter(Method == "CITS_spill_3pct")
p_trd_biasNHJ<-p_trd_biasNHJ%>%arrange(true_val,rho,n)%>%select(true_val,B_hat,Bias,rho,n,everything())

bs_bias_effect1 <- p_trd_biasNHJ %>%
  group_by(EstimandScenario) %>%
  summarise(
    mean_abs_bias =round(mean(Bias, na.rm = TRUE),2),
    mcse_bias=round(mean(MCSE_Bias, na.rm = TRUE),3),
    mean_percent_Bias=round(mean(Bias_percent,na.rm =TRUE),2),
    mean_estimate=round(mean(B_hat,na.rm=TRUE),4),
    Empirical_estimate=round(mean(Empirical_SE,na.rm=TRUE),3),
    .groups = "drop"
  )

bs_bias_effect1



bias_by_rho <- p_trd_biasNHJ %>%
  group_by(EstimandScenario,rho) %>%
  summarise(
    mean_abs_bias =round(mean(Bias, na.rm = TRUE),2),
    mcse_bias=round(mean(MCSE_Bias, na.rm = TRUE),3),
    mean_percent_Bias=round(mean(Bias_percent,na.rm =TRUE),2),
    mean_estimate=round(mean(B_hat,na.rm=TRUE),4),
    Empirical_estimateSE=round(mean(Empirical_SE,na.rm=TRUE),3),
    .groups = "drop"
  )

bias_by_rho
bias_by_rho$Empirical_estimateSE


bias_by_n <- p_trd_biasNHJ %>%
  group_by(EstimandScenario,n) %>%
  summarise(
    mean_abs_bias =round(mean(Bias, na.rm = TRUE),2),
    mcse_bias=round(mean(MCSE_Bias, na.rm = TRUE),3),
    mean_percent_Bias=round(mean(Bias_percent,na.rm =TRUE),2),
    mean_estimate=round(mean(B_hat,na.rm=TRUE),4),
    Empirical_estimateSE=round(mean(Empirical_SE,na.rm=TRUE),3),
    .groups = "drop"
  )

bias_by_n
bias_by_n$Empirical_estimateSE


bias_by_n %>%
  filter(
    EstimandScenario == "Moderate",
    n %in% c(12, 150)
  )


#Empirical vs modelbased SE
se_by_rho <- p_trd_biasNHJ %>%
  filter(EstimandScenario == "Moderate") %>%
  group_by(rho) %>%
  summarise(
    mean_empirical_SE = round(mean(Empirical_SE, na.rm = TRUE), 3),
    mean_model_SE = round(mean(avg_model_SE, na.rm = TRUE), 3),
    mean_SE_ratio = round(mean(ratio_Emp_ModelSE, na.rm = TRUE), 2),
    .groups = "drop"
  )

se_by_rho


se_by_rho %>%
  filter(rho %in% c(0, 0.4, 0.8))



se_by_n <- p_trd_biasNHJ %>%
  filter(
    EstimandScenario == "Moderate",
    rho %in% c(0.4, 0.8),
    n %in% c(12, 150)
  ) %>%
  select(
    rho,
    n,
    Empirical_SE,
    avg_model_SE,
    ratio_Emp_ModelSE
  ) %>%
  arrange(rho, n) %>%
  mutate(
    Empirical_SE = round(Empirical_SE, 3),
    avg_model_SE = round(avg_model_SE, 3),
    ratio_Emp_ModelSE = round(ratio_Emp_ModelSE, 2)
  )

se_by_n


se_by_effect_rho <- p_trd_biasNHJ %>%
  filter(rho %in% c(0, 0.4, 0.8)) %>%
  group_by(EstimandScenario, rho) %>%
  summarise(
    mean_empirical_SE = round(mean(Empirical_SE, na.rm = TRUE), 3),
    mean_model_SE = round(mean(avg_model_SE, na.rm = TRUE), 3),
    mean_SE_ratio = round(mean(ratio_Emp_ModelSE, na.rm = TRUE), 2),
    .groups = "drop"
  )

se_by_effect_rho




##Coverage

coverage_values <- p_trd_biasNHJ %>%
  filter(
    EstimandScenario == "Moderate",
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 80, 150)
  ) %>%
  select(
    rho,
    n,
    coverage_pct,
    MCSE_cov_pct
  ) %>%
  arrange(rho, n) %>%
  mutate(
    coverage_pct = round(coverage_pct, 2),
    MCSE_cov_pct = round(MCSE_cov_pct, 2)
  )

coverage_values

#POWER
power_by_effect <- p_trd_biasNHJ %>%
  group_by(EstimandScenario) %>%
  summarise(
    mean_power = round(mean(power_pct, na.rm = TRUE), 2),
    min_power  = round(min(power_pct, na.rm = TRUE), 2),
    max_power  = round(max(power_pct, na.rm = TRUE), 2),
    mean_MCSE  = round(mean(MCSE_power_pct, na.rm = TRUE), 3),
    .groups = "drop"
  )

power_by_effect


power_moderate <- p_trd_biasNHJ %>%
  filter(
    EstimandScenario == "Moderate",
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    rho,
    n,
    power_pct,
    MCSE_power_pct
  ) %>%
  arrange(rho, n) %>%
  mutate(
    power_pct = round(power_pct, 2),
    MCSE_power_pct = round(MCSE_power_pct, 3)
  )

power_moderate


power_threshold <- p_trd_biasNHJ %>%
  group_by(EstimandScenario, rho) %>%
  summarise(
    first_n_80 = ifelse(
      any(power_pct >= 80, na.rm = TRUE),
      min(n[power_pct >= 80], na.rm = TRUE),
      NA
    ),
    .groups = "drop"
  )

power_threshold


#mse
## =========================================================
## MSE summaries for quantitative interpretation
## =========================================================

# 1. Overall MSE by intervention effect size
mse_by_effect <- p_trd_biasNHJ %>%
  group_by(EstimandScenario) %>%
  summarise(
    mean_MSE = round(mean(MSE_estimate, na.rm = TRUE), 3),
    min_MSE  = round(min(MSE_estimate, na.rm = TRUE), 3),
    max_MSE  = round(max(MSE_estimate, na.rm = TRUE), 3),
    mean_MCSE_MSE = round(mean(MCSE_MSE, na.rm = TRUE), 3),
    .groups = "drop"
  )

mse_by_effect


# 2. Moderate effect: selected autocorrelation levels and series lengths
mse_moderate <- p_trd_biasNHJ %>%
  filter(
    EstimandScenario == "Moderate",
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    rho,
    n,
    MSE_estimate,
    MCSE_MSE
  ) %>%
  arrange(rho, n) %>%
  mutate(
    MSE_estimate = round(MSE_estimate, 3),
    MCSE_MSE = round(MCSE_MSE, 3)
  )

mse_moderate


# 3. Compare effect sizes at short and long series
mse_effect_compare <- p_trd_biasNHJ %>%
  filter(
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 150)
  ) %>%
  select(
    EstimandScenario,
    rho,
    n,
    MSE_estimate,
    MCSE_MSE
  ) %>%
  arrange(EstimandScenario, rho, n) %>%
  mutate(
    MSE_estimate = round(MSE_estimate, 3),
    MCSE_MSE = round(MCSE_MSE, 3)
  )

mse_effect_compare





##Objective 2:CITS vs multivariable regression method

## =========================================================
## Objective 2: Bias summaries for CITS vs Multivariable
## =========================================================

p_compare <- perf_plot %>%
  filter(Method %in% c("Trd", "CITS")) %>%
  arrange(Method, true_val, rho, n)


## 1. Mean estimate and mean bias by method and effect size
bias_ridge_compare <- p_compare %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    true_effect   = round(first(true_val), 4),
    mean_estimate = round(mean(B_hat, na.rm = TRUE), 4),
    mean_bias     = round(mean(Bias, na.rm = TRUE), 2),
    MCSE_bias     = round(mean(MCSE_Bias, na.rm = TRUE), 3),
    .groups = "drop"
  )


bias_ridge_compare %>%
  mutate(
    true_effect = sprintf("%.4f", true_effect),
    mean_estimate = sprintf("%.4f", mean_estimate),
    mean_bias = sprintf("%.2f", mean_bias),
    MCSE_bias = sprintf("%.3f", MCSE_bias)
  )
## 2. Bias by method, effect size and autocorrelation
## This checks whether the apparent centring changes as rho increases
bias_ridge_rho <- p_compare %>%
  group_by(Method, EstimandScenario, rho) %>%
  summarise(
    mean_estimate = round(mean(B_hat, na.rm = TRUE), 4),
    mean_bias     = round(mean(Bias, na.rm = TRUE), 2),
    MCSE_bias     = round(mean(MCSE_Bias, na.rm = TRUE), 3),
    .groups = "drop"
  )

bias_ridge_rho


##percent bias
## =========================================================
## OBJECTIVE 2: PERCENT BIAS
## Multivariable regression versus CITS
## =========================================================

p_compare_bias <- perf_plot %>%
  filter(Method %in% c("Trd", "CITS"))


## 1. Overall percent bias by method and intervention effect
percent_bias_effect <- p_compare_bias %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_percent_bias = round(mean(Bias_percent, na.rm = TRUE), 2),
    min_percent_bias  = round(min(Bias_percent, na.rm = TRUE), 2),
    max_percent_bias  = round(max(Bias_percent, na.rm = TRUE), 2),
    .groups = "drop"
  )

percent_bias_effect


## 2. Percent bias by method, effect size and autocorrelation
percent_bias_rho <- p_compare_bias %>%
  group_by(Method, EstimandScenario, rho) %>%
  summarise(
    mean_percent_bias = round(mean(Bias_percent, na.rm = TRUE), 2),
    .groups = "drop"
  )

percent_bias_rho


## 3. Selected values for the SMALL effect
## This is where the methods appear to differ most
percent_bias_small <- p_compare_bias %>%
  filter(
    EstimandScenario == "Small",
    rho %in% c(0, 0.4, 0.6, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    rho,
    n,
    Bias_percent
  ) %>%
  arrange(Method, rho, n) %>%
  mutate(
    Bias_percent = round(Bias_percent, 2)
  )

percent_bias_small


## 4. Moderate and large effects:
## compact overall comparison by method
percent_bias_modlarge <- p_compare_bias %>%
  filter(EstimandScenario %in% c("Moderate", "Large")) %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_percent_bias = round(mean(Bias_percent, na.rm = TRUE), 2),
    max_percent_bias  = round(max(Bias_percent, na.rm = TRUE), 2),
    .groups = "drop"
  )

percent_bias_modlarge





  ##SENSITIVITY
## =========================================================
## SENSITIVITY ANALYSIS
## Standard CITS versus control rho fixed at 0.4
## =========================================================

p_sensitivity <- perf_plot %>%
  filter(Method %in% c("CITS", "CITS_0.4_constant"))


## 1. Overall percent bias by method and effect size
sensitivity_effect <- p_sensitivity %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_percent_bias = round(mean(Bias_percent, na.rm = TRUE), 2),
    min_percent_bias  = round(min(Bias_percent, na.rm = TRUE), 2),
    max_percent_bias  = round(max(Bias_percent, na.rm = TRUE), 2),
    .groups = "drop"
  )

sensitivity_effect


## 2. Percent bias by outcome autocorrelation
sensitivity_rho <- p_sensitivity %>%
  group_by(Method, EstimandScenario, rho) %>%
  summarise(
    mean_percent_bias = round(mean(Bias_percent, na.rm = TRUE), 2),
    .groups = "drop"
  )

sensitivity_rho


## 3. Direct comparison at selected series lengths
sensitivity_selected <- p_sensitivity %>%
  filter(
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    EstimandScenario,
    rho,
    n,
    Bias_percent
  ) %>%
  arrange(EstimandScenario, rho, n, Method) %>%
  mutate(
    Bias_percent = round(Bias_percent, 2)
  )

sensitivity_selected



##Empirical vs model based standard error
 ##Empirical standard error

## =========================================================
## Objective 2: Empirical standard error
## CITS versus multivariable regression
## =========================================================

se_compare <- perf_plot %>%
  filter(Method %in% c("Trd", "CITS"))


## 1. Overall empirical SE by method and effect size
se_by_effect <- se_compare %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_emp_SE = round(mean(Empirical_SE, na.rm = TRUE), 3),
    min_emp_SE  = round(min(Empirical_SE, na.rm = TRUE), 3),
    max_emp_SE  = round(max(Empirical_SE, na.rm = TRUE), 3),
    .groups = "drop"
  )

se_by_effect


## 2. Moderate effect:
## empirical SE by autocorrelation
se_moderate_rho <- se_compare %>%
  filter(EstimandScenario == "Moderate") %>%
  group_by(Method, rho) %>%
  summarise(
    mean_emp_SE = round(mean(Empirical_SE, na.rm = TRUE), 3),
    .groups = "drop"
  )

se_moderate_rho


## 3. Moderate effect at selected series lengths
se_moderate_selected <- se_compare %>%
  filter(
    EstimandScenario == "Moderate",
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    rho,
    n,
    Empirical_SE
  ) %>%
  arrange(rho, n, Method) %>%
  mutate(
    Empirical_SE = round(Empirical_SE, 3)
  )

se_moderate_selected


## 4. Direct CITS / multivariable empirical-SE ratio
## Values < 1 indicate smaller empirical SE for CITS
se_ratio <- se_compare %>%
  filter(EstimandScenario == "Moderate") %>%
  select(Method, rho, n, Empirical_SE) %>%
  tidyr::pivot_wider(
    names_from = Method,
    values_from = Empirical_SE
  ) %>%
  mutate(
    CITS_to_MV_ratio = round(CITS / Trd, 2)
  ) %>%
  filter(
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  arrange(rho, n)

se_ratio




## =========================================================
## Objective 2:
## Empirical versus model-based SE
## Multivariable regression versus CITS
## =========================================================

se_alignment <- perf_plot %>%
  filter(Method %in% c("Trd", "CITS"))


## 1. Overall empirical and model-based SE by method/effect size
se_alignment_effect <- se_alignment %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_emp_SE   = round(mean(Empirical_SE, na.rm = TRUE), 3),
    mean_model_SE = round(mean(avg_model_SE, na.rm = TRUE), 3),
    mean_ratio    = round(mean(ratio_Emp_ModelSE, na.rm = TRUE), 2),
    .groups = "drop"
  )

se_alignment_effect


## 2. Moderate effect by autocorrelation
se_alignment_rho <- se_alignment %>%
  filter(EstimandScenario == "Moderate") %>%
  group_by(Method, rho) %>%
  summarise(
    mean_emp_SE   = round(mean(Empirical_SE, na.rm = TRUE), 3),
    mean_model_SE = round(mean(avg_model_SE, na.rm = TRUE), 3),
    mean_ratio    = round(mean(ratio_Emp_ModelSE, na.rm = TRUE), 2),
    .groups = "drop"
  )

se_alignment_rho


## 3. Moderate effect at selected series lengths
se_alignment_selected <- se_alignment %>%
  filter(
    EstimandScenario == "Moderate",
    rho %in% c(0, 0.4, 0.6, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    rho,
    n,
    Empirical_SE,
    avg_model_SE,
    ratio_Emp_ModelSE
  ) %>%
  arrange(Method, rho, n) %>%
  mutate(
    Empirical_SE      = round(Empirical_SE, 3),
    avg_model_SE      = round(avg_model_SE, 3),
    ratio_Emp_ModelSE = round(ratio_Emp_ModelSE, 2)
  )

se_alignment_selected






## Sensitivity analysis:
## Empirical SE for standard CITS versus CITS with control rho fixed at 0.4

se_sensitivity <- perf_plot %>%
  filter(Method %in% c("CITS", "CITS_0.4_constant"))

## 1. Moderate effect by outcome autocorrelation
se_sensitivity_rho <- se_sensitivity %>%
  filter(EstimandScenario == "Moderate") %>%
  group_by(Method, rho) %>%
  summarise(
    mean_emp_SE = round(mean(Empirical_SE, na.rm = TRUE), 3),
    .groups = "drop"
  )

se_sensitivity_rho


## 2. Moderate effect at selected series lengths
se_sensitivity_selected <- se_sensitivity %>%
  filter(
    EstimandScenario == "Moderate",
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    rho,
    n,
    Empirical_SE
  ) %>%
  arrange(rho, n, Method) %>%
  mutate(
    Empirical_SE = round(Empirical_SE, 3)
  )

se_sensitivity_selected




##Coverage

## =========================================================
## Coverage comparison: CITS vs Multivariable regression
## =========================================================

## =========================================================
## COVERAGE COMPARISON: CITS VS MULTIVARIABLE REGRESSION
## =========================================================

library(dplyr)

## Keep only the two main methods
coverage_dat <- perf_plot %>%
  filter(Method %in% c("Trd", "CITS"))


## ---------------------------------------------------------
## 1. Overall coverage by method and intervention effect size
## ---------------------------------------------------------

coverage_effect <- coverage_dat %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_coverage = round(mean(coverage_pct, na.rm = TRUE), 2),
    min_coverage  = round(min(coverage_pct, na.rm = TRUE), 2),
    max_coverage  = round(max(coverage_pct, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  as.data.frame()

coverage_effect


## ---------------------------------------------------------
## 2. Mean coverage by method, effect size and autocorrelation
## ---------------------------------------------------------

coverage_rho <- coverage_dat %>%
  group_by(Method, EstimandScenario, rho) %>%
  summarise(
    mean_coverage = round(mean(coverage_pct, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(EstimandScenario, rho, Method) %>%
  as.data.frame()

coverage_rho


## ---------------------------------------------------------
## 3. Coverage at selected series lengths
## ---------------------------------------------------------

coverage_selected <- coverage_dat %>%
  filter(
    rho %in% c(0, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    EstimandScenario,
    rho,
    n,
    coverage_pct
  ) %>%
  arrange(
    EstimandScenario,
    rho,
    n,
    Method
  ) %>%
  mutate(
    coverage_pct = round(coverage_pct, 2)
  ) %>%
  as.data.frame()

coverage_selected




##power
## =========================================================
## POWER COMPARISON: CITS VS MULTIVARIABLE REGRESSION
## =========================================================

library(dplyr)

## Keep only the two main methods
power_dat <- perf_plot %>%
  filter(Method %in% c("Trd", "CITS"))


## ---------------------------------------------------------
## 1. Overall power by method and intervention effect size
## ---------------------------------------------------------

power_effect <- power_dat %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_power = round(mean(power_pct, na.rm = TRUE), 2),
    min_power  = round(min(power_pct, na.rm = TRUE), 2),
    max_power  = round(max(power_pct, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  as.data.frame()

power_effect


## ---------------------------------------------------------
## 2. Mean power by method, effect size and autocorrelation
## ---------------------------------------------------------

power_rho <- power_dat %>%
  group_by(Method, EstimandScenario, rho) %>%
  summarise(
    mean_power = round(mean(power_pct, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(EstimandScenario, rho, Method) %>%
  as.data.frame()

power_rho


## ---------------------------------------------------------
## 3. Power at selected series lengths
## ---------------------------------------------------------

power_selected <- power_dat %>%
  filter(
    rho %in% c(0, 0.2, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    EstimandScenario,
    rho,
    n,
    power_pct
  ) %>%
  arrange(
    EstimandScenario,
    rho,
    n,
    Method
  ) %>%
  mutate(
    power_pct = round(power_pct, 2)
  ) %>%
  as.data.frame()

power_selected



## =========================================================
## MSE COMPARISON: CITS VS MULTIVARIABLE REGRESSION
## =========================================================

library(dplyr)

## Keep only the two main methods
mse_dat <- perf_plot %>%
  filter(Method %in% c("Trd", "CITS"))


## ---------------------------------------------------------
## 1. Overall MSE by method and intervention effect size
## ---------------------------------------------------------

mse_effect <- mse_dat %>%
  group_by(Method, EstimandScenario) %>%
  summarise(
    mean_MSE = round(mean(MSE_estimate, na.rm = TRUE), 3),
    min_MSE  = round(min(MSE_estimate, na.rm = TRUE), 3),
    max_MSE  = round(max(MSE_estimate, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  as.data.frame()

mse_effect


## ---------------------------------------------------------
## 2. Mean MSE by method, effect size and autocorrelation
## ---------------------------------------------------------

mse_rho <- mse_dat %>%
  group_by(Method, EstimandScenario, rho) %>%
  summarise(
    mean_MSE = round(mean(MSE_estimate, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(EstimandScenario, rho, Method) %>%
  as.data.frame()

mse_rho


## ---------------------------------------------------------
## 3. MSE at selected series lengths
## ---------------------------------------------------------

mse_selected <- mse_dat %>%
  filter(
    rho %in% c(0, 0.2, 0.4, 0.8),
    n %in% c(12, 40, 80, 150)
  ) %>%
  select(
    Method,
    EstimandScenario,
    rho,
    n,
    MSE_estimate
  ) %>%
  arrange(
    EstimandScenario,
    rho,
    n,
    Method
  ) %>%
  mutate(
    MSE_estimate = round(MSE_estimate, 3)
  ) %>%
  as.data.frame()

mse_selected












##Abstract extract
##uNDER CONTAMINATION
p_trd_biasNHJ<- perf_plot %>%filter(Method == "CITS_spill_3pct")
summary_bias <- p_trd_biasNHJ %>%
  group_by(EstimandScenario) %>%
  summarise(
    mean_percent_bias = round(mean(Bias_percent, na.rm = TRUE),2),
    sd_percent_bias   = round(sd(Bias_percent, na.rm = TRUE),2)
  )

summary_bias

p_trd_biasNHJ1 <- perf_plot %>% 
  filter(Method == "CITS_spill_3pct") %>% 
  mutate(
    autocorrelation_group = case_when(
      rho <= 0.4 ~ "Low-to-moderate (rho <= 0.4)",
      rho == 0.8 ~ "High (rho = 0.8)",
      TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(autocorrelation_group))

summary_performance <- p_trd_biasNHJ1 %>% 
  group_by(
    EstimandScenario,
    autocorrelation_group
  ) %>% 
  summarise(
    mean_percent_bias = round(
      mean(Bias_percent, na.rm = TRUE), 2
    ),
    
    sd_percent_bias = round(
      sd(Bias_percent, na.rm = TRUE), 2
    ),
    
    mean_coverage = round(
      mean(coverage_pct, na.rm = TRUE), 2
    ),
    
    mean_empirical_SE = round(
      mean(Empirical_SE, na.rm = TRUE), 3
    ),
    
    .groups = "drop"
  )

summary_performance



#cits VS MULTIVARIABLE
perf_plot %>%
  filter(
    Method %in% c("Trd", "CITS"),
    EstimandScenario == "Moderate",
    rho %in% c(0, 0.2, 0.4, 0.8)
  ) %>%
  mutate(
    autocorrelation = case_when(
      rho <= 0.4 ~ "Low-to-moderate autocorrelation",
      rho == 0.8 ~ "Strong autocorrelation"
    )
  ) %>%
  group_by(Method, autocorrelation) %>%
  summarise(
    mean_bias = mean(Bias_percent, na.rm = TRUE),
    sd_bias = sd(Bias_percent, na.rm = TRUE),
    min_bias = min(Bias_percent, na.rm = TRUE),
    max_bias = max(Bias_percent, na.rm = TRUE),
    min_EmpSE = min(Empirical_SE, na.rm = TRUE),
    max_EmpSE = max(Empirical_SE, na.rm = TRUE),
    min_coverage = min(coverage_pct, na.rm = TRUE),
    max_coverage = max(coverage_pct, na.rm = TRUE),
    .groups = "drop"
  )
