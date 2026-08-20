library(mice)
library(dplyr)
library(gtsummary)
library(sandwich)
library(MASS)
library(knitr)
##Read data for the imputation of birth weight and analysis
FINAL_DATA_08_Apr_2026_u<-readRDS("data/FINAL_DATA_08_Apr_2026_u.rds")
FINAL_DATA_08_Apr_2026<-FINAL_DATA_08_Apr_2026_u
tab <- table(FINAL_DATA_08_Apr_2026$birth_weight, useNA = "ifany")
prop.table(tab)
summary(FINAL_DATA_08_Apr_2026$birth_weight)

nrow(FINAL_DATA_08_Apr_2026[is.na(FINAL_DATA_08_Apr_2026$birth_weight),])
#prop missing
print(round((nrow(FINAL_DATA_08_Apr_2026[is.na(FINAL_DATA_08_Apr_2026$birth_weight), ]) /
               nrow(FINAL_DATA_08_Apr_2026)) * 100, 2))


# 1. Prepare imputation dataset
imp_data <- FINAL_DATA_08_Apr_2026 %>%
  mutate(
    # clean impossible birth weight codes
    #birth_weight = ifelse(birth_weight >= 9999, NA, birth_weight),,
    
    # make categorical variables factors
    #new_Location_birth = as.factor(new_Location_birth),
    #Child_month_birthFactor = as.factor(Child_month_birthFactor),
    
    # keep year as continuous/numeric
    child_year_birth = as.numeric(child_year_birth),
    dist_kchKm = as.numeric(dist_kchKm),
    maternal_age_yrs = as.numeric(maternal_age_yrs),
    birth_weightorign=as.numeric(birth_weight)
  ) %>%
  dplyr::select(
    birth_weight,
    birth_weightorign
    ,
    new_Location_birth,
    dist_kchKm,
    child_year_birth,
    Child_month_birthFactor,
    maternal_age_yrs
  )


imp_data <- imp_data %>%
  mutate(
    new_Location_birth = as.character(new_Location_birth),
    new_Location_birth = ifelse(new_Location_birth == "Missing",
                                NA,
                                new_Location_birth),
    new_Location_birth = factor(new_Location_birth)
  )

##lets do some cleaning of impossible values before imputation
imp_data <- imp_data %>%
  mutate(
    birth_weight = case_when(
      birth_weight < 500 ~ NA_real_,
      birth_weight > 6000 ~ NA_real_,
      TRUE ~ birth_weight
    )
  )

#prop missing
print(round((nrow(imp_data[is.na(imp_data$birth_weight), ]) /
               nrow(imp_data)) * 100, 2))



ini <- mice(imp_data, maxit = 0)

meth <- ini$method
pred <- ini$predictorMatrix

# Start by not imputing anything
meth[] <- ""

# Impute variables needed so birth_weight can be fully imputed
meth["birth_weight"] <- "pmm"
meth["new_Location_birth"] <- "polyreg"
meth["maternal_age_yrs"] <- "pmm"

# Predictor matrix
pred[,] <- 0

# Impute birth_weight using all predictors
pred["birth_weight", c(
  "new_Location_birth",
  "dist_kchKm",
  "child_year_birth",
  "Child_month_birthFactor",
  "maternal_age_yrs"
)] <- 1

# Impute new_Location_birth using other variables
pred["new_Location_birth", c(
  "birth_weight",
  "dist_kchKm",
  "child_year_birth",
  "Child_month_birthFactor",
  "maternal_age_yrs"
)] <- 1

# Impute maternal_age_yrs using other variables
pred["maternal_age_yrs", c(
  "birth_weight",
  "new_Location_birth",
  "dist_kchKm",
  "child_year_birth",
  "Child_month_birthFactor"
)] <- 1

set.seed(12345)

imp_birthweight <- mice(
  imp_data,
  m = 20,
  method = meth,
  predictorMatrix = pred,
  maxit = 20,
  donors = 5,
  printFlag = TRUE
)

imp_long <- complete(
  imp_birthweight,
  action = "long",
  include = FALSE
)

imp_long$dataset_id <- paste0("dataset_", imp_long$.imp)

##check missing birth weight
View(imp_long %>%
  group_by(.imp) %>%
  summarise(
    n_total = n(),
    n_missing_birth_weight = sum(is.na(birth_weight)),
    percent_missing_birth_weight =
      round(mean(is.na(birth_weight)) * 100, 2),
    mean_birth_weight = round(mean(birth_weight, na.rm = TRUE), 2),
    median_birth_weight = round(median(birth_weight, na.rm = TRUE), 2)
  ))

View(imp_long)

imp_long$LBW<-ifelse(imp_long$birth_weight<2500,1,0)
imp_long$LBW1<-ifelse(imp_long$birth_weightorign<2500,1,0)
imp_longPLay<-imp_long
imp_longPLay<-imp_longPLay[!is.na(imp_longPLay$birth_weightorign),]
##lets do some summary ststistics
bw_summary_by_imp <- imp_long %>%
  group_by(.imp, dataset_id) %>%
  summarise(
    mean_birth_weight = round(mean(birth_weight, na.rm = TRUE), 2),
    median_birth_weight = round(median(birth_weight, na.rm = TRUE), 2),
    min_birth_weight = min(birth_weight, na.rm = TRUE),
    max_birth_weight = max(birth_weight, na.rm = TRUE),
    
    mean_birth_weight_orig = round(mean(birth_weightorign, na.rm = TRUE), 2),
    median_birth_weight_orig = round(median(birth_weightorign, na.rm = TRUE), 2),
    min_birth_weight_orig = min(birth_weightorign, na.rm = TRUE),
    max_birth_weight_orig = max(birth_weightorign, na.rm = TRUE),
    
    .groups = "drop"
  )

#bw_summary_by_imp

lbw_monthly <- imp_long %>%
  group_by(.imp, dataset_id, child_year_birth, Child_month_birthFactor) %>%
  summarise(
    total_births = n(),
    number_lbw = sum(LBW == 1, na.rm = TRUE),
    incidence_lbw_per_100000 =
      round((number_lbw / total_births) * 100000, 2),
    .groups = "drop"
  ) %>%
  rename(
    yoa = child_year_birth,
    admission_month = Child_month_birthFactor
  ) %>%
  mutate(
    admission_month = match(
      tolower(as.character(admission_month)),
      tolower(month.name)
    )
  ) %>%
  arrange(.imp, yoa, admission_month)

View(lbw_monthly)

lbw_monthly<-lbw_monthly[!(lbw_monthly$yoa==2011 & lbw_monthly$admission_month<4),]
lbw_monthly<-lbw_monthly[lbw_monthly$yoa<=2019,]

saveRDS(lbw_monthly,"data/lbw_monthlyr.rds")


##Lets read admissions
AdmissionInc<-readRDS("data/AdmissionInc.rds")
AdmissionIncUpdated_Final<-left_join(AdmissionInc,lbw_monthly,by =c("yoa","admission_month"))
AdmissionIncUpdated_Final$total_births<-NULL
AdmissionIncUpdated_Final$incidence_lbw_per_1000<-(AdmissionIncUpdated_Final$number_lbw/AdmissionIncUpdated_Final$LiveBirths_currentdata)*1000
AdmissionIncUpdated_Final<-AdmissionIncUpdated_Final%>%arrange(.imp,yoa,admission_month)
AdmissionIncUpdated_Final<-AdmissionIncUpdated_Final %>%
  dplyr::select(.imp, yoa, admission_month,LiveBirths,LiveBirths_currentdata,number_lbw, everything())

##lets prepare our data for analysis
##generate time variable
AdmissionIncUpdated_Final <- AdmissionIncUpdated_Final %>%
  group_by(.imp) %>%
  arrange(yoa, admission_month, .by_group = TRUE) %>%
  mutate(
    time = row_number(),
    Policy_Indicator = if_else(
      yoa > 2013 | (yoa == 2013 & admission_month >= 7),
      1, 0
    ),
    Policy_Indicator = if_else(
      yoa == 2013 & admission_month == 6,
      NA_real_,
      as.numeric(Policy_Indicator)
    ),
    Post_Interraction=(time-28)*Policy_Indicator
  ) %>%
  ungroup()


AdmissionIncUpdated_Final<-AdmissionIncUpdated_Final%>%mutate(Health_WorkerStrike=ifelse((AdmissionIncUpdated_Final$yoa==2011 & admission_month==12)|(AdmissionIncUpdated_Final$yoa==2012 & admission_month %in% c(3,9,10,12))|(AdmissionIncUpdated_Final$yoa==2013 & admission_month %in% c(1,2,12))|(AdmissionIncUpdated_Final$yoa==2016 & admission_month==12),1,0))

AdmissionIncUpdated_Final<-AdmissionIncUpdated_Final %>%
  dplyr::select(.imp, yoa, admission_month,time,Policy_Indicator,Post_Interraction,Health_WorkerStrike,LiveBirths,LiveBirths_currentdata,number_lbw, everything())


View(AdmissionIncUpdated_Final)


AdmissionIncUpdated_Final<-AdmissionIncUpdated_Final%>%
  dplyr::arrange(time) %>%
  dplyr::mutate(
    month_dummy= factor(admission_month, levels = 1:12, labels = month.abb),
    #t_c   = t - mean(t, na.rm = TRUE),
    #off   = log(control + delta)     # <-- offset with small positive delta
  )

AdmissionIncUpdated_Final<-AdmissionIncUpdated_Final%>%arrange(.imp,yoa,admission_month)


##checking if control is nice(My control real seem to be working.See unaffected by policy and meeting common trend assumption)
## Pre-policy data (You can use any of the imputed datasets)
##common trend assumption test
df_pre <- AdmissionIncUpdated_Final %>%
  filter(Policy_Indicator == 0, .imp ==2)

m_nb <- MASS::glm.nb(
  deathsunder1 ~ time + month_dummy + offset(log(deaths5to14)),
  data = df_pre
)

n <- nrow(df_pre)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["time"]
se_time <- sqrt(diag(V_hac))["time"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

common_trend_table <- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 2)
  ) %>%
  dplyr::select(IRR, `95% CI`, `p-value`)

common_trend_table
##Is my control affected by the policy;You Can any of th imputed dataset
##check both the immediate impact and change in trend
dt<-AdmissionIncUpdated_Final[AdmissionIncUpdated_Final$.imp==1,]

n <- nrow(dt)
L <- floor(n^(1/4))


m_ctrl <- MASS::glm.nb(
  deaths5to14 ~ time+Post_Interraction+ Policy_Indicator+ month_dummy+offset(log(PersonMonth5to14_All)),
  data =dt
)
summary(m_ctrl)


## Coefficients of interest


# IMPORTANT: Recalculate HAC covariance matrix
V_hac <- sandwich::NeweyWest(
  m_ctrl,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)
terms <- c("Policy_Indicator", "Post_Interraction")

## Extract coefficients and HAC standard errors
beta <- coef(m_ctrl)[terms]
se   <- sqrt(diag(V_hac))[terms]

## Wald statistics and p-values
z <- beta / se
p <- 2 * (1 - pnorm(abs(z)))

## IRRs and 95% confidence intervals
results_control <- data.frame(
  Effect = c(
    "Immediate level change",
    "Change in post-intervention trend"
  ),
  IRR = exp(beta),
  lower = exp(beta - 1.96 * se),
  upper = exp(beta + 1.96 * se),
  p_value = p
)

## Format results
results_control <- results_control %>%
  mutate(
    IRR = sprintf("%.2f", IRR),
    `95% CI` = paste0(
      sprintf("%.2f", lower),
      "–",
      sprintf("%.2f", upper)
    ),
    `p-value` = sprintf("%.3f", p_value)
  ) %>%
  dplyr::select(
    Effect,
    IRR,
    `95% CI`,
    `p-value`
  )

results_control




   ##creating data to test DID
   #dtanaDID<-dtana
   #dtanaDID$.imp<-NULL
    #dtanaDID$dataset_id<-NULL
   ##treated
    #dtanaDIDtreat<-dtanaDID%>%dplyr::select(yoa,admission_month,month_dummy,Health_WorkerStrike,deathsunder1)%>%rename(count=deathsunder1)
    #dtanaDIDtreat<-dtanaDIDtreat%>%mutate(treat=as.integer(!is.na(yoa)))
    #control
    #dtanaDIDcontrol<-dtanaDID%>%dplyr::select(yoa,admission_month,month_dummy,Health_WorkerStrike,deaths5to14)%>%rename(count=deaths5to14)
    #dtanaDIDcontrol<-dtanaDIDcontrol%>%mutate(treat=0)
    #Append datasets
    #dtanaDID1<-rbind(dtanaDIDtreat,dtanaDIDcontrol)
    #dtanaDID1<-dtanaDID1%>%mutate(Policy_Indicator = if_else(
      #yoa > 2013 | (yoa == 2013 & admission_month >= 7),
      #1, 0
    #))
    
    #dtanaDID1$Policy_Indicator<-ifelse(dtanaDID1$yoa==2013 & dtanaDID1$admission_month==6,NA_real_,as.numeric(dtanaDID1$Policy_Indicator))
     #building DID model
    #m_nb <- MASS::glm.nb(
      #count ~ treat * Policy_Indicator,
      #data =dtanaDID1
    #)
    
    ##lets addd healtcare woker strike and monthly dummy in DID nad see
    #m_nbdid<- MASS::glm.nb(
      #count ~ treat * Policy_Indicator+month_dummy+Health_WorkerStrike,
      #data =dtanaDID1
    #)
    
    
    
    #ndid<- nrow(dtanaDID1)
    #ndid
    #Ldid<- floor(n^(1/4))
    
    
    
    ##Doing the analysis now. Using control series:Post_Interraction
    ##Linda mama was from April 2017
AdmissionIncUpdated_Final<-AdmissionIncUpdated_Final%>%mutate(LindaMama=case_when(yoa>2017|(yoa==2017 & admission_month>=4)~1,TRUE~0))

    dtana<-AdmissionIncUpdated_Final[AdmissionIncUpdated_Final$.imp==3,]
    
    write.csv(
      dtana,
      "data/dtana.csv",
      row.names = FALSE
    )
n <- nrow(dtana)
n
L <- floor(n^(1/4))

##Descriptive Table
library(dplyr)
library(tidyr)

#===========================================================
# 1. DEFINE PRE- AND POST-INTERVENTION PERIODS
#===========================================================

dtana_desc <- dtana %>%
  mutate(
    Period = case_when(
      
      # Pre-intervention: April 2011 to May 2013
      (yoa == 2011 & admission_month >= 4) |
        (yoa == 2012) |
        (yoa == 2013 & admission_month <= 5) ~
        "Pre-intervention",
      
      # Post-intervention: July 2013 to December 2019
      (yoa == 2013 & admission_month >= 7) |
        (yoa >= 2014 & yoa <= 2019) ~
        "Post-intervention",
      
      # June 2013 excluded
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period)) %>%
  mutate(
    Period = factor(
      Period,
      levels = c("Pre-intervention", "Post-intervention")
    )
  )


#===========================================================
# 2. CALCULATE MONTHLY MORTALITY RATE PER 1,000 PERSON-YEARS
#===========================================================

dtana_desc <- dtana_desc %>%
  mutate(
    mortality_rate_1000 =
      (deathsunder1 /person_Month_under1) * 1000
  )


#===========================================================
# 3. CALCULATE MEAN, SD, MEDIAN AND IQR BY PERIOD
#===========================================================

summary_stats <- dtana_desc %>%
  group_by(Period) %>%
  summarise(
    
    # Death counts
    death_mean = mean(deathsunder1, na.rm = TRUE),
    death_sd = sd(deathsunder1, na.rm = TRUE),
    death_median = median(deathsunder1, na.rm = TRUE),
    death_q1 = quantile(deathsunder1, 0.25, na.rm = TRUE),
    death_q3 = quantile(deathsunder1, 0.75, na.rm = TRUE),
    
    # Mortality rates
    rate_mean = mean(mortality_rate_1000, na.rm = TRUE),
    rate_sd = sd(mortality_rate_1000, na.rm = TRUE),
    rate_median = median(mortality_rate_1000, na.rm = TRUE),
    rate_q1 = quantile(mortality_rate_1000, 0.25, na.rm = TRUE),
    rate_q3 = quantile(mortality_rate_1000, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )


#===========================================================
# 4. CREATE EXACT TABLE STRUCTURE
#===========================================================

table_final <- data.frame(
  
  Characteristic = c(
    "Mean (SD)",
    "Median (IQR)"
  ),
  
  `Case counts: Pre-intervention` = c(
    sprintf(
      "%.2f (%.2f)",
      summary_stats$death_mean[
        summary_stats$Period == "Pre-intervention"
      ],
      summary_stats$death_sd[
        summary_stats$Period == "Pre-intervention"
      ]
    ),
    
    sprintf(
      "%.2f (%.2f–%.2f)",
      summary_stats$death_median[
        summary_stats$Period == "Pre-intervention"
      ],
      summary_stats$death_q1[
        summary_stats$Period == "Pre-intervention"
      ],
      summary_stats$death_q3[
        summary_stats$Period == "Pre-intervention"
      ]
    )
  ),
  
  `Case counts: Post-intervention` = c(
    sprintf(
      "%.2f (%.2f)",
      summary_stats$death_mean[
        summary_stats$Period == "Post-intervention"
      ],
      summary_stats$death_sd[
        summary_stats$Period == "Post-intervention"
      ]
    ),
    
    sprintf(
      "%.2f (%.2f–%.2f)",
      summary_stats$death_median[
        summary_stats$Period == "Post-intervention"
      ],
      summary_stats$death_q1[
        summary_stats$Period == "Post-intervention"
      ],
      summary_stats$death_q3[
        summary_stats$Period == "Post-intervention"
      ]
    )
  ),
  
  `Mortality rate per 1,000: Pre-intervention` = c(
    sprintf(
      "%.2f (%.2f)",
      summary_stats$rate_mean[
        summary_stats$Period == "Pre-intervention"
      ],
      summary_stats$rate_sd[
        summary_stats$Period == "Pre-intervention"
      ]
    ),
    
    sprintf(
      "%.2f (%.2f–%.2f)",
      summary_stats$rate_median[
        summary_stats$Period == "Pre-intervention"
      ],
      summary_stats$rate_q1[
        summary_stats$Period == "Pre-intervention"
      ],
      summary_stats$rate_q3[
        summary_stats$Period == "Pre-intervention"
      ]
    )
  ),
  
  `Mortality rate per 1,000: Post-intervention` = c(
    sprintf(
      "%.2f (%.2f)",
      summary_stats$rate_mean[
        summary_stats$Period == "Post-intervention"
      ],
      summary_stats$rate_sd[
        summary_stats$Period == "Post-intervention"
      ]
    ),
    
    sprintf(
      "%.2f (%.2f–%.2f)",
      summary_stats$rate_median[
        summary_stats$Period == "Post-intervention"
      ],
      summary_stats$rate_q1[
        summary_stats$Period == "Post-intervention"
      ],
      summary_stats$rate_q3[
        summary_stats$Period == "Post-intervention"
      ]
    )
  ),
  
  check.names = FALSE
)

table_final


library(dplyr)
library(gt)
library(flextable)
library(officer)
library(knitr)
library(kableExtra)

table_gt <- table_final %>%
  gt() %>%
  
  # ----------------------------------------------------------
# Rename the individual columns
# ----------------------------------------------------------
cols_label(
  `Characteristic` = "Characteristic",
  `Case counts: Pre-intervention` = "Pre-intervention",
  `Case counts: Post-intervention` = "Post-intervention",
  `Mortality rate per 1,000: Pre-intervention` = "Pre-intervention",
  `Mortality rate per 1,000: Post-intervention` = "Post-intervention"
) %>%
  
  # ----------------------------------------------------------
# Group the columns exactly like your handwritten table
# ----------------------------------------------------------
tab_spanner(
  label = "Case counts",
  columns = c(
    `Case counts: Pre-intervention`,
    `Case counts: Post-intervention`
  )
) %>%
  
  tab_spanner(
    label = "Mortality rate per 1,000 person-years",
    columns = c(
      `Mortality rate per 1,000: Pre-intervention`,
      `Mortality rate per 1,000: Post-intervention`
    )
  ) %>%
  
  # ----------------------------------------------------------
# Table title
# ----------------------------------------------------------
tab_header(
  title = md("**Under-1 mortality before and after policy implementation**"),
  subtitle = md(
    "Monthly case counts and crude mortality rates, April 2011–December 2019"
  )
) %>%
  
  # ----------------------------------------------------------
# Footnote
# ----------------------------------------------------------
tab_source_note(
  source_note = md(
    "Pre-intervention: April 2011–May 2013; 
       post-intervention: July 2013–December 2019. 
       June 2013 was excluded as the intervention month. 
       Mortality rates are expressed per 1,000 under-1 person-years."
  )
) %>%
  
  # ----------------------------------------------------------
# Alignment
# ----------------------------------------------------------
cols_align(
  align = "left",
  columns = Characteristic
) %>%
  
  cols_align(
    align = "center",
    columns = -Characteristic
  ) %>%
  
  # ----------------------------------------------------------
# Formatting
# ----------------------------------------------------------
tab_style(
  style = cell_text(weight = "bold"),
  locations = cells_column_spanners()
) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(11),
    data_row.padding = px(6),
    column_labels.padding = px(6),
    table.width = pct(100)
  )

# ============================================================
# VIEW TABLE
# ============================================================

table_gt


##latex code
latex_table <- table_final

# ============================================================
# Generate LaTeX
# ============================================================

latex_code <- latex_table %>%
  
  kbl(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", "c", "c", "c", "c"),
    
    # IMPORTANT:
    # Use only the label name here.
    # Do NOT write "tab:under1_descriptive",
    # otherwise you can get tab:tab:under1_descriptive.
    label = "under1_descriptive",
    
    caption = paste0(
      "\\textbf{Under-1 mortality before and after policy implementation}",
      "\\\\",
      "\\normalfont Monthly case counts and crude mortality rates, ",
      "April 2011--December 2019"
    ),
    
    col.names = c(
      "\\textbf{Characteristic}",
      "\\makecell[c]{\\textbf{Pre-}\\\\\\textbf{intervention}}",
      "\\makecell[c]{\\textbf{Post-}\\\\\\textbf{intervention}}",
      "\\makecell[c]{\\textbf{Pre-}\\\\\\textbf{intervention}}",
      "\\makecell[c]{\\textbf{Post-}\\\\\\textbf{intervention}}"
    )
  ) %>%
  
  # ==========================================================
# GROUPED HEADERS
# ==========================================================

add_header_above(
  c(
    " " = 1,
    "\\textbf{Case counts}" = 2,
    "\\makecell[c]{\\textbf{Mortality rate per 1,000 person-}\\\\\\textbf{years}}" = 2
  ),
  escape = FALSE,
  line = TRUE
) %>%
  
  # ==========================================================
# WIDTH / POSITION FORMATTING
# ==========================================================

kable_styling(
  latex_options = c(
    "hold_position"
  ),
  full_width = FALSE,
  position = "center"
) %>%
  
  # Characteristic column
  column_spec(
    1,
    width = "2.7cm"
  ) %>%
  
  # Case count columns
  column_spec(
    2,
    width = "3.0cm"
  ) %>%
  
  column_spec(
    3,
    width = "3.0cm"
  ) %>%
  
  # Mortality-rate columns
  column_spec(
    4,
    width = "3.0cm"
  ) %>%
  
  column_spec(
    5,
    width = "3.0cm"
  ) %>%
  
  # ==========================================================
# TABLE NOTE
# ==========================================================

footnote(
  general = paste0(
    "Pre-intervention: April 2011--May 2013; ",
    "post-intervention: July 2013--December 2019. ",
    "June 2013 was excluded as the intervention month. ",
    "Mortality rates are expressed per 1,000 under-1 person-years."
  ),
  general_title = "",
  threeparttable = TRUE,
  escape = FALSE
)


# ============================================================
# VIEW THE GENERATED LATEX CODE
# ============================================================

latex_code



writeLines(
  as.character(latex_code),
  "data/Table1_under1_mortality.tex"
)



## Boxplot
dtana_desc <- dtana_desc %>%
  mutate(
    Period_plot = factor(
      Period,
      levels = c("Pre-intervention", "Post-intervention")
    )
  )

# Create boxplot
p_deaths_box <- ggplot(
  dtana_desc,
  aes(
    x = Period_plot,
    y = deathsunder1,
    fill = Period_plot
  )
) +
  
  geom_boxplot(
    width = 0.55,
    linewidth = 0.8,
    staplewidth = 0.35,
    outlier.shape = 21,
    outlier.size = 2.8,
    outlier.stroke = 0.8
  ) +
  
  # Blue = pre-intervention; red = post-intervention
  scale_fill_manual(
    values = c(
      "Pre-intervention" = "#6BAED6",
      "Post-intervention" = "#FB6A4A"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Monthly under-1 deaths"
  ) +
  
  guides(fill = "none") +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    panel.grid = element_blank()
  )

# View plot
p_deaths_box


# Save as PDF
ggsave(
  filename = "graphs/Under1_deaths_boxplot.pdf",
  plot = p_deaths_box,
  width = 7,
  height = 5,
  units = "in"
)

# Save as high-resolution PNG
ggsave(
  filename = "graphs/Under1_deaths_boxplot.png",
  plot = p_deaths_box,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600
)



##lets try full model with interraction fast
m_nb <- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+Post_Interraction+month_dummy+Health_WorkerStrike+offset(log(deaths5to14)),
  data =dtana)
summary(m_nb)

m_nb <- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+month_dummy+Health_WorkerStrike+ offset(log(deaths5to14)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Policy_Indicator"]
se_time <- sqrt(diag(V_hac))["Policy_Indicator"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 3)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS




##Now lets build the multivariable regression model men?
##I will rely on the first imputed dataset to decide on the predictors then build the same model across all imputations I dont want to completicate MSc work as if am earning nobel price.
dtana<-AdmissionIncUpdated_Final[AdmissionIncUpdated_Final$.imp==1,]
#dtana$MortalityPer1000<-(dtana$deathsunder1/dtana$person_Month_under1)*1000

   ##Start with Univariate model
##Linda mama programme  PersonMonthunder1_All
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+LindaMama+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["LindaMama"]
se_time <- sqrt(diag(V_hac))["LindaMama"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 2)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS




#m_nb<- MASS::glm.nb(
  #deathsunder1 ~Policy_Indicator+LindaMama+ offset(log(PersonMonthunder1_All)),
  #data =dtana
#)

#n <- nrow(dtana)
#L <- floor(n^(1/4))

#V_hac <- sandwich::NeweyWest(
  #m_nb,
  #lag = L,
  #prewhite = FALSE,
 # adjust = TRUE
#)

#beta_time <- coef(m_nb)["LindaMama"]
#se_time <- sqrt(diag(V_hac))["LindaMama"]

#z_time <- beta_time / se_time
#p_time <- 2 * (1 - pnorm(abs(z_time)))

#policy_effect_tableCITS<- data.frame(
  #IRR = exp(beta_time),
  #lower = exp(beta_time - 1.96 * se_time),
  #upper = exp(beta_time + 1.96 * se_time),
  #p_value = p_time
#) %>%
  #mutate(
    #IRR = round(IRR, 2),
    #`95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    #`p-value` = round(p_value, 2)
  #) %>%dplyr::select(IRR, `95% CI`, `p-value`)

#policy_effect_tableCITS






    #Health_WorkerStrike 0.02;
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+Health_WorkerStrike+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Health_WorkerStrike"]
se_time <- sqrt(diag(V_hac))["Health_WorkerStrike"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 2)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS
  





#measles(Exclude)
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+Measles_Incidence+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Measles_Incidence"]
se_time <- sqrt(diag(V_hac))["Measles_Incidence"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 2)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS

#HIV(0.0071);Not sure but may include it in the multivariable model
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+HIV_Incidence+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["HIV_Incidence"]
se_time <- sqrt(diag(V_hac))["HIV_Incidence"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 2)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS

#Malariaincidence (Include)
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+Malaria_Incidence+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Malaria_Incidence"]
se_time <- sqrt(diag(V_hac))["Malaria_Incidence"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 4),
    `95% CI` = paste0(round(lower, 4), "–", round(upper, 4)),
    `p-value` = round(p_value, 4)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS


#neonatal sepsis incidence (include p=0.014)
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+neonatal_sepsis_Incidence+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["neonatal_sepsis_Incidence"]
se_time <- sqrt(diag(V_hac))["neonatal_sepsis_Incidence"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 4),
    `95% CI` = paste0(round(lower, 4), "–", round(upper, 4)),
    `p-value` = round(p_value, 4)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS



##Congenital abnormalitie (Include p-value=0.0943)
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+birth_defects_Incidence+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["birth_defects_Incidence"]
se_time <- sqrt(diag(V_hac))["birth_defects_Incidence"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 4),
    `95% CI` = paste0(round(lower, 4), "–", round(upper, 4)),
    `p-value` = round(p_value, 4)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS


##Birth Asphia (Exclude p-value=0.5886)
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+birth_Asphyxia_Incidence+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["birth_Asphyxia_Incidence"]
se_time <- sqrt(diag(V_hac))["birth_Asphyxia_Incidence"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 4),
    `95% CI` = paste0(round(lower, 4), "–", round(upper, 4)),
    `p-value` = round(p_value, 4)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS



#Malnutration (Exclude p.value=0.1213 ) only including variable from univariate with p-value <0.1
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+malnutrition_Incidence+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["malnutrition_Incidence"]
se_time <- sqrt(diag(V_hac))["malnutrition_Incidence"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 4),
    `95% CI` = paste0(round(lower, 4), "–", round(upper, 4)),
    `p-value` = round(p_value, 4)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS



#Low birth weiht incidence (EXCLUDE P=0.2146);eCLUDING THIS VARIABLE IS EVEN NICE FOR ME AS i WILL NOT USE THE IMPUTED DATASETS HENCE AVoiding RUBIN combination
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+incidence_lbw_per_100000+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["incidence_lbw_per_100000"]
se_time <- sqrt(diag(V_hac))["incidence_lbw_per_100000"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 4),
    `95% CI` = paste0(round(lower, 4), "–", round(upper, 4)),
    `p-value` = round(p_value, 4)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS


# Time/yoa admission (Inlude p=0.0619)
m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+yoa+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["yoa"]
se_time <- sqrt(diag(V_hac))["yoa"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 4),
    `95% CI` = paste0(round(lower, 4), "–", round(upper, 4)),
    `p-value` = round(p_value, 4)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS

##Variables for multivariable before using likelihood ratio test
#Health_WorkerStrike+HIV_Incidence+Malaria_Incidence+neonatal_sepsis_Incidence+birth_defects_Incidence+yoa+month_dummy;Am tempted to try inluding LBW incidence;INCIDENCE OF LOW BIRTH WAIT MAKE irr=1.203 exclude it kabisa


#Checking overdispersion
library(performance)
##lets fit a poisson to asses the overdispersion
# Poisson multivariable model

# Overdispersion test

#dispersion ratio =   2.265
#Pearson's Chi-Squared = 192.549
                #p-value = < 0.001

#Overdispersion detected.

m_pois <- glm(
  deathsunder1 ~ Policy_Indicator +
    Health_WorkerStrike +
    HIV_Incidence +
    Malaria_Incidence +
    neonatal_sepsis_Incidence +
    birth_defects_Incidence +
    yoa +
    month_dummy +
    offset(log(person_Month_under1)),
  family = poisson(link = "log"),
  data = dtana
)

check_overdispersion(m_pois)



##fit negative binomial as we have overdispersion

m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+Health_WorkerStrike+HIV_Incidence+Malaria_Incidence+neonatal_sepsis_Incidence+birth_defects_Incidence+yoa+month_dummy+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Policy_Indicator"]
se_time <- sqrt(diag(V_hac))["Policy_Indicator"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 3)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS


summary(m_nb)
##lets build full model for appendix



library(dplyr)
library(stringr)
library(sandwich)
library(knitr)
library(kableExtra)

#--------------------------------------------------
# 1. Newey-West HAC variance-covariance matrix
#--------------------------------------------------

n <- nobs(m_nb)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

#--------------------------------------------------
# 2. Extract coefficients and HAC standard errors
#--------------------------------------------------

beta <- coef(m_nb)
se_hac <- sqrt(diag(V_hac))

z_hac <- beta / se_hac
p_hac <- 2 * pnorm(-abs(z_hac))

#--------------------------------------------------
# 3. Calculate IRRs and 95% CIs
#--------------------------------------------------

full_multivariable <- data.frame(
  term = names(beta),
  IRR = exp(beta),
  lower = exp(beta - 1.96 * se_hac),
  upper = exp(beta + 1.96 * se_hac),
  p_value = p_hac,
  row.names = NULL
)

#--------------------------------------------------
# 4. Give variables reader-friendly names
#--------------------------------------------------

full_multivariable <- full_multivariable %>%
  mutate(
    Variable = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "Policy_Indicator" ~ "Free Maternity Policy",
      term == "Health_WorkerStrike" ~ "Healthcare-worker strike period",
      term == "HIV_Incidence" ~ "HIV incidence",
      term == "Malaria_Incidence" ~ "Malaria incidence",
      term == "neonatal_sepsis_Incidence" ~ "Neonatal sepsis incidence",
      term == "birth_defects_Incidence" ~ "Birth defects incidence",
      term == "yoa" ~ "Year of admission",
      
      term == "month_dummyFeb" ~ "February",
      term == "month_dummyMar" ~ "March",
      term == "month_dummyApr" ~ "April",
      term == "month_dummyMay" ~ "May",
      term == "month_dummyJun" ~ "June",
      term == "month_dummyJul" ~ "July",
      term == "month_dummyAug" ~ "August",
      term == "month_dummySep" ~ "September",
      term == "month_dummyOct" ~ "October",
      term == "month_dummyNov" ~ "November",
      term == "month_dummyDec" ~ "December",
      
      TRUE ~ term
    )
  )

#--------------------------------------------------
# 5. Add January as reference category
#--------------------------------------------------

jan_reference <- data.frame(
  term = "month_dummyJan",
  IRR = 1,
  lower = NA,
  upper = NA,
  p_value = NA,
  Variable = "January (Reference)"
)

# Place January immediately before February
month_position <- which(full_multivariable$term == "month_dummyFeb")

full_multivariable <- bind_rows(
  full_multivariable[1:(month_position - 1), ],
  jan_reference,
  full_multivariable[month_position:nrow(full_multivariable), ]
)

#--------------------------------------------------
# 6. Format for thesis
#--------------------------------------------------

full_multivariable_table <- full_multivariable %>%
  mutate(
    IRR = ifelse(
      Variable == "January (Reference)",
      "1.00",
      sprintf("%.2f", IRR)
    ),
    
    `95% CI` = ifelse(
      Variable == "January (Reference)",
      "Reference",
      paste0(
        sprintf("%.2f", lower),
        "--",
        sprintf("%.2f", upper)
      )
    ),
    
    `p-value` = case_when(
      Variable == "January (Reference)" ~ "---",
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) %>%dplyr::select(
    Variable,
    IRR,
    `95% CI`,
    `p-value`
  )

full_multivariable_table



latex_multivariable <- full_multivariable_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", "c", "c", "c"),
    caption = "Full results from the multivariable negative binomial regression model",
    label = "full_multivariable"
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )

latex_multivariable
   ##lets test interraction
m_nbInte<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+Health_WorkerStrike+HIV_Incidence+Malaria_Incidence+neonatal_sepsis_Incidence+birth_defects_Incidence+yoa+month_dummy+ offset(log(person_Month_under1)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nbInte,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Policy_Indicator"]
se_time <- sqrt(diag(V_hac))["Policy_Indicator"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 3)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS















##CITS model
#fitPoisson model and cehck overdispersion
m_pois <- glm(
  deathsunder1 ~ Policy_Indicator +
    Health_WorkerStrike +month_dummy +
    offset(log(deaths5to14)),
  family = poisson(link = "log"),
  data = dtana
)

m_nb<-m_pois

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Policy_Indicator"]
se_time <- sqrt(diag(V_hac))["Policy_Indicator"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 3)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS

# CITS model:Overdispesion detected:P-value
#dispersion ratio =   5.378
#Pearson's Chi-Squared = 483.991
                #p-value = < 0.001

#Overdispersion detected.

check_overdispersion(m_pois)




##Fit negative binomial model as there was an evidence of overdispersion
m_nb <- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+month_dummy+Health_WorkerStrike+ offset(log(deaths5to14)),
  data =dtana
)

n <- nrow(dtana)
L <- floor(n^(1/4))

V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

beta_time <- coef(m_nb)["Policy_Indicator"]
se_time <- sqrt(diag(V_hac))["Policy_Indicator"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

policy_effect_tableCITS<- data.frame(
  IRR = exp(beta_time),
  lower = exp(beta_time - 1.96 * se_time),
  upper = exp(beta_time + 1.96 * se_time),
  p_value = p_time
) %>%
  mutate(
    IRR = round(IRR, 2),
    `95% CI` = paste0(round(lower, 2), "–", round(upper, 2)),
    `p-value` = round(p_value, 3)
  ) %>%dplyr::select(IRR, `95% CI`, `p-value`)

policy_effect_tableCITS

m_cits<-m_nb
m_cits




#--------------------------------------------------
# 1. Fit CITS negative binomial model full model
#--------------------------------------------------

m_cits <- MASS::glm.nb(
  deathsunder1 ~ Policy_Indicator +
    month_dummy +
    Health_WorkerStrike +
    offset(log(deaths5to14)),
  data = dtana
)

summary(m_cits)

#--------------------------------------------------
# 2. Newey-West HAC variance-covariance matrix
#--------------------------------------------------

n <- nobs(m_cits)
L <- floor(n^(1/4))

V_hac_cits <- sandwich::NeweyWest(
  m_cits,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

#--------------------------------------------------
# 3. Extract coefficients and HAC standard errors
#--------------------------------------------------

beta_cits <- coef(m_cits)
se_hac_cits <- sqrt(diag(V_hac_cits))

z_hac_cits <- beta_cits / se_hac_cits
p_hac_cits <- 2 * pnorm(-abs(z_hac_cits))

#--------------------------------------------------
# 4. Calculate IRRs and 95% confidence intervals
#--------------------------------------------------

full_cits <- data.frame(
  term = names(beta_cits),
  IRR = exp(beta_cits),
  lower = exp(beta_cits - 1.96 * se_hac_cits),
  upper = exp(beta_cits + 1.96 * se_hac_cits),
  p_value = p_hac_cits,
  row.names = NULL
)

#--------------------------------------------------
# 5. Give variables reader-friendly names
#--------------------------------------------------

full_cits <- full_cits %>%
  mutate(
    Variable = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "Policy_Indicator" ~ "Free Maternity Policy",
      term == "Health_WorkerStrike" ~ "Healthcare-worker strike period",
      
      term == "month_dummyFeb" ~ "February",
      term == "month_dummyMar" ~ "March",
      term == "month_dummyApr" ~ "April",
      term == "month_dummyMay" ~ "May",
      term == "month_dummyJun" ~ "June",
      term == "month_dummyJul" ~ "July",
      term == "month_dummyAug" ~ "August",
      term == "month_dummySep" ~ "September",
      term == "month_dummyOct" ~ "October",
      term == "month_dummyNov" ~ "November",
      term == "month_dummyDec" ~ "December",
      
      TRUE ~ term
    )
  )

#--------------------------------------------------
# 6. Add January as reference category
#--------------------------------------------------

jan_reference_cits <- data.frame(
  term = "month_dummyJan",
  IRR = 1,
  lower = NA,
  upper = NA,
  p_value = NA,
  Variable = "January (Reference)"
)

month_position <- which(full_cits$term == "month_dummyFeb")

full_cits <- bind_rows(
  full_cits[1:(month_position - 1), ],
  jan_reference_cits,
  full_cits[month_position:nrow(full_cits), ]
)

#--------------------------------------------------
# 7. Format for thesis
#--------------------------------------------------

full_cits_table <- full_cits %>%
  mutate(
    IRR = ifelse(
      Variable == "January (Reference)",
      "1.00",
      sprintf("%.2f", IRR)
    ),
    
    `95% CI` = ifelse(
      Variable == "January (Reference)",
      "Reference",
      paste0(
        sprintf("%.2f", lower),
        "--",
        sprintf("%.2f", upper)
      )
    ),
    
    `p-value` = case_when(
      Variable == "January (Reference)" ~ "---",
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) %>%
  dplyr::select(
    Variable,
    IRR,
    `95% CI`,
    `p-value`
  )

full_cits_table

#--------------------------------------------------
# 8. Create LaTeX-ready table
#--------------------------------------------------

latex_cits <- full_cits_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", "c", "c", "c"),
    caption = "Full results from the controlled interrupted time series negative binomial model",
    label = "full_cits"
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )

latex_cits













#--------------------------------------------------
# 1. Multivariable NB model WITH interaction
#--------------------------------------------------

m_nb_int <- MASS::glm.nb(
  deathsunder1 ~ Policy_Indicator +
    Post_Interraction +
    Health_WorkerStrike +
    HIV_Incidence +
    Malaria_Incidence +
    neonatal_sepsis_Incidence +
    birth_defects_Incidence +
    yoa +
    month_dummy +
    offset(log(person_Month_under1)),
  data = dtana
)

n_nb <- nobs(m_nb_int)
L_nb <- floor(n_nb^(1/4))

V_hac_nb <- sandwich::NeweyWest(
  m_nb_int,
  lag = L_nb,
  prewhite = FALSE,
  adjust = TRUE
)

beta_nb <- coef(m_nb_int)["Post_Interraction"]
se_nb <- sqrt(diag(V_hac_nb))["Post_Interraction"]

z_nb <- beta_nb / se_nb
p_nb <- 2 * pnorm(-abs(z_nb))

#--------------------------------------------------
# 2. CITS NB model WITH interaction
#--------------------------------------------------

m_cits_int <- MASS::glm.nb(
  deathsunder1 ~ Policy_Indicator +
    Post_Interraction +
    month_dummy +
    Health_WorkerStrike +
    offset(log(deaths5to14)),
  data = dtana
)

n_cits <- nobs(m_cits_int)
L_cits <- floor(n_cits^(1/4))

V_hac_cits <- sandwich::NeweyWest(
  m_cits_int,
  lag = L_cits,
  prewhite = FALSE,
  adjust = TRUE
)

beta_cits <- coef(m_cits_int)["Post_Interraction"]
se_cits <- sqrt(diag(V_hac_cits))["Post_Interraction"]

z_cits <- beta_cits / se_cits
p_cits <- 2 * pnorm(-abs(z_cits))

#--------------------------------------------------
# 3. Combine interaction results
#--------------------------------------------------

interaction_table <- data.frame(
  Model = c(
    "Multivariable regression",
    "Controlled interrupted time series"
  ),
  IRR = c(
    exp(beta_nb),
    exp(beta_cits)
  ),
  lower = c(
    exp(beta_nb - 1.96 * se_nb),
    exp(beta_cits - 1.96 * se_cits)
  ),
  upper = c(
    exp(beta_nb + 1.96 * se_nb),
    exp(beta_cits + 1.96 * se_cits)
  ),
  p_value = c(
    p_nb,
    p_cits
  )
) %>%
  mutate(
    IRR = sprintf("%.2f", IRR),
    `95% CI` = paste0(
      sprintf("%.2f", lower),
      "--",
      sprintf("%.2f", upper)
    ),
    `p-value` = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  ) %>%
  dplyr::select(
    Model,
    IRR,
    `95% CI`,
    `p-value`
  )

interaction_table


latex_interaction <- interaction_table %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = c("l", "c", "c", "c"),
    caption = "Assessment of post-intervention trend change in the empirical models",
    label = "interaction_empirical"
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )

latex_interaction









#Model diagonistics

m_cits <- MASS::glm.nb(
  deathsunder1 ~ Policy_Indicator +
    month_dummy +
    Health_WorkerStrike +
    offset(log(deaths5to14)),
  data = dtana
)

m_nb<- MASS::glm.nb(
  deathsunder1 ~Policy_Indicator+Health_WorkerStrike+HIV_Incidence+Malaria_Incidence+neonatal_sepsis_Incidence+birth_defects_Incidence+yoa+month_dummy+ offset(log(person_Month_under1)),
  data =dtana
)


library(ggplot2)
library(gridExtra)

#--------------------------------------------------
# 1. Extract Pearson residuals
#--------------------------------------------------

res_mult <- residuals(m_nb, type = "pearson")
res_cits <- residuals(m_cits, type = "pearson")

#--------------------------------------------------
# 2. ACF and PACF objects
#--------------------------------------------------

acf_mult <- acf(
  res_mult,
  lag.max = 12,
  plot = FALSE
)

pacf_mult <- pacf(
  res_mult,
  lag.max = 12,
  plot = FALSE
)

acf_cits <- acf(
  res_cits,
  lag.max = 12,
  plot = FALSE
)

pacf_cits <- pacf(
  res_cits,
  lag.max = 12,
  plot = FALSE
)

#--------------------------------------------------
# 3. Function to convert ACF/PACF to ggplot
#--------------------------------------------------

make_acf_plot <- function(acf_obj, title, n_obs, ylab = "ACF") {
  
  df <- data.frame(
    lag = as.numeric(acf_obj$lag),
    value = as.numeric(acf_obj$acf)
  )
  
  ggplot(df, aes(x = lag, y = value)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend = lag, yend = 0)) +
    geom_hline(
      yintercept = c(
        1.96 / sqrt(n_obs),
        -1.96 / sqrt(n_obs)
      ),
      linetype = "dashed"
    ) +
    labs(
      title = title,
      x = "Lag (months)",
      y = ylab
    ) +
    scale_x_continuous(breaks = 0:12) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}


make_pacf_plot <- function(pacf_obj, title, n_obs) {
  
  df <- data.frame(
    lag = as.numeric(pacf_obj$lag),
    value = as.numeric(pacf_obj$acf)
  )
  
  ggplot(df, aes(x = lag, y = value)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend = lag, yend = 0)) +
    geom_hline(
      yintercept = c(
        1.96 / sqrt(n_obs),
        -1.96 / sqrt(n_obs)
      ),
      linetype = "dashed"
    ) +
    labs(
      title = title,
      x = "Lag (months)",
      y = "PACF"
    ) +
    scale_x_continuous(breaks = 1:12) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

#--------------------------------------------------
# 4. Create the four plots
#--------------------------------------------------

n_mult <- nobs(m_nb)
n_cits <- nobs(m_cits)

p1 <- make_acf_plot(
  acf_mult,
  #A) Multivariable regression: ACF",
  "A",
  n_mult
)

p2 <- make_pacf_plot(
  pacf_mult,
  #B) Multivariable regression: PACF",
  "B",
  n_mult
)

p3 <- make_acf_plot(
  acf_cits,
  #"C) CITS: ACF",
  "C",
  n_cits
)

p4 <- make_pacf_plot(
  pacf_cits,
  #"D) CITS: PACF",
  "D",
  n_cits
)

autocorrelation_plot <- grid.arrange(
  p1, p2, p3, p4,
  ncol = 2
)

ggsave(
  "graphs/autocorrelation_diagnostics.png",
  autocorrelation_plot,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  "graphs/autocorrelation_diagnostics.pdf",
  autocorrelation_plot,
  width = 10,
  height = 8
)



##Ljung test

#--------------------------------------------------
# Ljung-Box test for residual autocorrelation
#--------------------------------------------------

lb_mult <- Box.test(
  residuals(m_nb, type = "pearson"),
  lag = 12,
  type = "Ljung-Box"
)

lb_cits <- Box.test(
  residuals(m_cits, type = "pearson"),
  lag = 12,
  type = "Ljung-Box"
)

lb_mult
lb_cits


#==================================================
# AUTOMATICALLY GENERATE AUTOCORRELATION TABLE
#==================================================

# Ljung-Box tests
lb_mult <- Box.test(
  residuals(m_nb, type = "pearson"),
  lag = 12,
  type = "Ljung-Box"
)

lb_cits <- Box.test(
  residuals(m_cits, type = "pearson"),
  lag = 12,
  type = "Ljung-Box"
)

#--------------------------------------------------
# Function for formatting p-values
#--------------------------------------------------

format_p <- function(p) {
  if (p < 0.001) {
    return("< 0.001")
  } else {
    return(sprintf("%.3f", p))
  }
}

#--------------------------------------------------
# Create table data
#--------------------------------------------------

autocorrelation_table <- data.frame(
  Model = c(
    "Multivariable regression",
    "Controlled interrupted time series"
  ),
  
  `Ljung--Box Q` = c(
    unname(lb_mult$statistic),
    unname(lb_cits$statistic)
  ),
  
  `Degrees of freedom` = c(
    unname(lb_mult$parameter),
    unname(lb_cits$parameter)
  ),
  
  `p-value` = c(
    format_p(lb_mult$p.value),
    format_p(lb_cits$p.value)
  ),
  
  check.names = FALSE
)

#--------------------------------------------------
# Round test statistics
#--------------------------------------------------

autocorrelation_table$`Ljung--Box Q` <-
  sprintf(
    "%.2f",
    as.numeric(autocorrelation_table$`Ljung--Box Q`)
  )

#--------------------------------------------------
# Create LaTeX rows
#--------------------------------------------------

latex_rows <- apply(
  autocorrelation_table,
  1,
  function(x) {
    paste(
      x[1],
      x[2],
      x[3],
      x[4],
      sep = " & "
    )
  }
)

latex_rows <- paste0(
  latex_rows,
  " \\\\"
)

#--------------------------------------------------
# Create complete LaTeX table
#--------------------------------------------------

latex_table <- paste0(
  "\\begin{table}[!h]
\\centering
\\caption{\\label{tab:autocorrelation}
\\textbf{Assessment of residual autocorrelation in the empirical models}}
\\begin{tabular}{lrrr}
\\toprule
Model & Ljung--Box $Q$ & Degrees of freedom & $p$-value \\\\
\\midrule
",
  paste(latex_rows, collapse = "\n"),
  "
\\bottomrule
\\end{tabular}

\\vspace{0.1cm}

\\begin{minipage}{0.95\\textwidth}
\\small
\\textit{Note:} The Ljung--Box test was applied to Pearson residuals
through lag 12. The null hypothesis was that the residual
autocorrelations through lag 12 were jointly equal to zero.
A $p$-value greater than 0.05 indicates no statistically significant
evidence of residual autocorrelation.
\\end{minipage}

\\end{table}
"
)

#--------------------------------------------------
# Write LaTeX file
#--------------------------------------------------

writeLines(
  latex_table,
  "data/autocorrelation_testljung.tex"
)

cat(latex_table)











