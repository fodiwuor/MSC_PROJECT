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
df_pre <- AdmissionIncUpdated_Final %>%
  filter(Policy_Indicator == 0, .imp == 2)

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
dt<-AdmissionIncUpdated_Final[AdmissionIncUpdated_Final$.imp==2,]

n <- nrow(dt)
L <- floor(n^(1/4))


m_ctrl <- MASS::glm.nb(
  deaths5to14 ~ time+ Policy_Indicator+ month_dummy,
  data =dt
)


## Newey–West variance–covariance matrix
V_hac <- sandwich::NeweyWest(
  m_ctrl,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)


## ---- Bottomley-style Wald test ----
beta_time <- coef(m_ctrl)["Policy_Indicator"]
se_time   <- sqrt(diag(V_hac))["Policy_Indicator"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

## Results
beta_time
se_time
z_time
p_time






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

    dtana<-AdmissionIncUpdated_Final[AdmissionIncUpdated_Final$.imp==2,]   
    
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

##CITS model
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

