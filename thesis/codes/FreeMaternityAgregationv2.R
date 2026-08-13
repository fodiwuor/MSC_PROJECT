##Read my admissions data
require(haven)
library(tidyverse)
library(lubridate)
library(stringr)
library(readxl)
library(MASS)       # glm.nb
library(sandwich)   # NeweyWest
library(lmtest)     # coeftest
require(dplyr)
AmissionMcsProjectFredMainAnalysis_StrikePeriodsRemoved <- read_dta("data/AmissionMcsProjectFredMainAnalysis_StrikePeriodsRemoved.dta")
admissionData<-AmissionMcsProjectFredMainAnalysis_StrikePeriodsRemoved
admissionData$admission_month <- month(admissionData$doa)
##sort based on date
admissionData<-admissionData%>%arrange(doa)
admissionData$Agecat<-ifelse(admissionData$agem<12,1,0)
#Agecat
##keep infants
admissionData<-admissionData[admissionData$Agecat==1,]
View(admissionData)
nrow(admissionData)

View(admissionData %>%dplyr::arrange(Agecat)%>%dplyr::select(serial,fk_person,doa,agem,agey,Agecat, everything()))
admissionData<-admissionData%>%dplyr::select(-Agecat)
View(admissionData)


table(admissionData$resident,useNA ="ifany")
##Incidence of malaria admissions:Aggregate malaria counts
Mal <- admissionData %>%
  dplyr::filter(resident == 1) %>%
  dplyr::mutate(mps2 = ifelse(is.na(mps), 99, mps)) %>%   # 99 = Missing
  dplyr::count(mps2, yoa, admission_month, name = "n") %>%
  complete(
    mps2 = c(0, 1, 9, 99),
    yoa = 2011:2019,
    admission_month = 1:12,
    fill = list(n = 0)
  )

View(
  admissionData[
    admissionData$yoa == 2011 &
      admissionData$admission_month == 4 &
      admissionData$resident == 1,
    c("yoa", "admission_month", "mps")
  ]
)

Mal<-Mal%>%dplyr::filter(mps2==1)%>%dplyr::rename(Malaria_count=n)%>%dplyr::select(-mps2)
Mal<-Mal%>%dplyr::filter(yoa>=2011 & yoa <=2019)
Mal<-Mal%>%dplyr::filter(!(yoa == 2011 & admission_month < 4))

##HIV related admissions
##hivresult
table(admissionData$hivresult,useNA ="ifany")
HIV<- admissionData %>%
  dplyr::filter(resident == 1) %>%
  dplyr::mutate(hivresult2= ifelse(is.na(hivresult), 99, hivresult)) %>%   # 99 = Missing
  dplyr::count(hivresult2, yoa, admission_month, name = "n") %>%
  complete(
    hivresult2 = c(0, 1, 9, 99),
    yoa = 2011:2019,
    admission_month = 1:12,
    fill = list(n = 0)
  )

View(
  admissionData[
    admissionData$yoa == 2011 &
      admissionData$admission_month == 4 &
      admissionData$resident == 1,
    c("yoa", "admission_month", "hivresult")
  ]
)

HIV<-HIV%>%dplyr::filter(hivresult2==1)%>%dplyr::rename(HIV_count=n)%>%dplyr::select(-hivresult2)
HIV<-HIV%>%dplyr::filter(yoa>=2011 & yoa <=2019)
HIV<-HIV%>%dplyr::filter(!(yoa == 2011 & admission_month < 4))
##count if HIV related admission equal expected totals;Nice it matches
HIV%>% summarise(total = sum(HIV_count)) %>% pull(total)


##Measles admissions
table(admissionData$Measles,useNA ="ifany")
Measles<- admissionData %>%
  dplyr::filter(resident == 1) %>%
  dplyr::mutate(Measles2= ifelse(is.na(Measles), 99, Measles)) %>%   # 99 = Missing
  dplyr::count(Measles2, yoa, admission_month, name = "n") %>%
  complete(
    Measles2 = c(0, 1,99),
    yoa = 2011:2019,
    admission_month = 1:12,
    fill = list(n = 0)
  )

View(
  admissionData[
    admissionData$yoa == 2011 &
      admissionData$admission_month == 4 &
      admissionData$resident == 1,
    c("yoa", "admission_month", "Measles")
  ]
)

Measles<-Measles%>%dplyr::filter(Measles2==1)%>%dplyr::rename(Measles_count=n)%>%dplyr::select(-Measles2)
Measles<-Measles%>%dplyr::filter(yoa>=2011 & yoa <=2019)
Measles<-Measles%>%dplyr::filter(!(yoa == 2011 & admission_month < 4))
##count if Measles related admission equal expected totals;Nice it matches
Measles%>% summarise(total = sum(Measles_count)) %>% pull(total)

##Neonatal sepsis admissions
table(admissionData$neonatal_sepsi,useNA ="ifany")
neonatal_sepsi<- admissionData %>%
  dplyr::filter(resident == 1) %>%
  dplyr::mutate(neonatal_sepsi2= ifelse(is.na(neonatal_sepsi), 99, neonatal_sepsi)) %>%   # 99 = Missing
  dplyr::count(neonatal_sepsi2, yoa, admission_month, name = "n") %>%
  complete(
    neonatal_sepsi2 = c(0, 1,99),
    yoa = 2011:2019,
    admission_month = 1:12,
    fill = list(n = 0)
  )

View(
  admissionData[
    admissionData$yoa == 2011 &
      admissionData$admission_month == 4 &
      admissionData$resident == 1,
    c("yoa", "admission_month", "neonatal_sepsi")
  ]
)

neonatal_sepsi<-neonatal_sepsi%>%dplyr::filter(neonatal_sepsi2==1)%>%dplyr::rename(neonatal_sepsis_count=n)%>%dplyr::select(-neonatal_sepsi2)
neonatal_sepsi<-neonatal_sepsi%>%dplyr::filter(yoa>=2011 & yoa <=2019)
neonatal_sepsi<-neonatal_sepsi%>%dplyr::filter(!(yoa == 2011 & admission_month < 4))
##count if Measles related admission equal expected totals;Nice it matches
neonatal_sepsi%>% summarise(total = sum(neonatal_sepsis_count)) %>% pull(total)


##congenital abnormalities/birth_defects
#name of variable is birth_defects
table(admissionData$birth_defects,useNA ="ifany")
birth_defects<- admissionData %>%
  dplyr::filter(resident == 1) %>%
  dplyr::mutate(birth_defects2= ifelse(is.na(birth_defects), 99, birth_defects)) %>%   # 99 = Missing
  dplyr::count(birth_defects2, yoa, admission_month, name = "n") %>%
  complete(
    birth_defects2 = c(0, 1,99),
    yoa = 2011:2019,
    admission_month = 1:12,
    fill = list(n = 0)
  )

View(
  admissionData[
    admissionData$yoa == 2011 &
      admissionData$admission_month == 4 &
      admissionData$resident == 1,
    c("yoa", "admission_month", "birth_defects")
  ]
)

birth_defects<-birth_defects%>%dplyr::filter(birth_defects2==1)%>%dplyr::rename(birth_defects_count=n)%>%dplyr::select(-birth_defects2)
birth_defects<-birth_defects%>%dplyr::filter(yoa>=2011 & yoa <=2019)
birth_defects<-birth_defects%>%dplyr::filter(!(yoa == 2011 & admission_month < 4))
##count if Measles related admission equal expected totals;Nice it matches
birth_defects%>% summarise(total = sum(birth_defects_count)) %>% pull(total)


##birth Asphyxia
table(admissionData$birth_asphyxia,useNA ="ifany")
birth_Asphyxia<- admissionData %>%
  dplyr::filter(resident == 1) %>%
  dplyr::mutate(birth_asphyxia2= ifelse(is.na(birth_asphyxia), 99, birth_asphyxia)) %>%   # 99 = Missing
  dplyr::count(birth_asphyxia2, yoa, admission_month, name = "n") %>%
  complete(
    birth_asphyxia2= c(0, 1,99),
    yoa = 2011:2019,
    admission_month = 1:12,
    fill = list(n = 0)
  )

View(
  admissionData[
    admissionData$yoa == 2011 &
      admissionData$admission_month == 4 &
      admissionData$resident == 1,
    c("yoa", "admission_month", "birth_asphyxia")
  ]
)

birth_Asphyxia<-birth_Asphyxia%>%dplyr::filter(birth_asphyxia2==1)%>%dplyr::rename(birth_Asphyxia_count=n)%>%dplyr::select(-birth_asphyxia2)
birth_Asphyxia<-birth_Asphyxia%>%dplyr::filter(yoa>=2011 & yoa <=2019)
birth_Asphyxia<-birth_Asphyxia%>%dplyr::filter(!(yoa == 2011 & admission_month < 4))
##count if birth_Asphyxia related admission equal expected totals;Nice it matches
birth_Asphyxia%>% summarise(total = sum(birth_Asphyxia_count)) %>% pull(total)

##Malnutrition admissions
table(admissionData$msmn,useNA ="ifany")
msmn<- admissionData %>%
  dplyr::filter(resident == 1) %>%
  dplyr::mutate(msmn2= ifelse(is.na(msmn), 99, msmn)) %>%   # 99 = Missing
  dplyr::count(msmn2, yoa, admission_month, name = "n") %>%
  complete(
    msmn2= c(0, 1,99),
    yoa = 2011:2019,
    admission_month = 1:12,
    fill = list(n = 0)
  )

View(
  admissionData[
    admissionData$yoa == 2011 &
      admissionData$admission_month == 4 &
      admissionData$resident == 1,
    c("yoa", "admission_month", "msmn")
  ]
)

msmn<-msmn%>%dplyr::filter(msmn2==1)%>%dplyr::rename(malnutrition_count=n)%>%dplyr::select(-msmn2)
msmn<-msmn%>%dplyr::filter(yoa>=2011 & yoa <=2019)
msmn<-msmn%>%dplyr::filter(!(yoa == 2011 & admission_month< 4))
##count if birth_Asphyxia related admission equal expected totals;Nice it matches
msmn%>% summarise(total = sum(malnutrition_count)) %>% pull(total)

##Merge the admissions data;
AdmissionInc<-inner_join(Mal, HIV, by =c("yoa","admission_month"))
AdmissionInc<-inner_join(AdmissionInc, Measles, by =c("yoa","admission_month"))
AdmissionInc<-inner_join(AdmissionInc, neonatal_sepsi, by =c("yoa","admission_month"))
AdmissionInc<-inner_join(AdmissionInc, birth_defects, by =c("yoa","admission_month"))
AdmissionInc<-inner_join(AdmissionInc, birth_Asphyxia, by =c("yoa","admission_month"))
AdmissionInc<-inner_join(AdmissionInc, msmn, by =c("yoa","admission_month"))
view(AdmissionInc)
head(AdmissionInc,n=10)

##Now lets read mid_year population table
temp14frematunder1 <- read_dta("data/temp14frematunder1.dta")
AdmissionInc<-inner_join(AdmissionInc, temp14frematunder1, by =c("yoa"))
AdmissionInc$person_Month<-((AdmissionInc$midyrpopunder1)/12)
AdmissionInc$person_Month5to14<-((AdmissionInc$midyrpop5to14)/12)
  ##adjust each month independently (Try this as best approach and check your estimate)
calendar_days_2011_2019Cht <- read.csv("C:/Users/THINKPAD/Desktop/DesktopDamp/MSC_STATISTICSNOTES/Personal_lerning notes/MSC_Project/documents/referencesTo_ReadForMethodFavouring/data/calendar_days_2011_2019Cht.csv")

calendar_days_2011_2019Cht<-calendar_days_2011_2019Cht%>%mutate(DayObserved=days_in_month) 

calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2011 &
    calendar_days_2011_2019Cht$month == "December"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2011 &
      calendar_days_2011_2019Cht$month == "December"
  ] - 9


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2012 &
    calendar_days_2011_2019Cht$month == "March"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2012 &
      calendar_days_2011_2019Cht$month == "March"
  ] - 15


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2012 &
    calendar_days_2011_2019Cht$month == "September"
] <-12


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2012 &
    calendar_days_2011_2019Cht$month == "October"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2012 &
      calendar_days_2011_2019Cht$month == "October"
  ] -4



calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2012 &
    calendar_days_2011_2019Cht$month == "December"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2012 &
      calendar_days_2011_2019Cht$month == "December"
  ] -22


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2013 &
    calendar_days_2011_2019Cht$month == "January"
] <-15

calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2013 &
    calendar_days_2011_2019Cht$month == "February"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2013 &
      calendar_days_2011_2019Cht$month == "February"
  ] -11


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2013 &
    calendar_days_2011_2019Cht$month == "December"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2013 &
      calendar_days_2011_2019Cht$month == "December"
  ] -14



calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2016 &
    calendar_days_2011_2019Cht$month == "December"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2016 &
      calendar_days_2011_2019Cht$month == "December"
  ] -27


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "January"
] <-0



calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "February"
] <-0





calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "March"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2017 &
      calendar_days_2011_2019Cht$month == "March"
  ] -15


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "June"
] <-4

calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "July"
] <-0

calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "August"
] <-0


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "September"
] <-0

calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "October"
] <-0


calendar_days_2011_2019Cht$DayObserved[
  calendar_days_2011_2019Cht$year == 2017 &
    calendar_days_2011_2019Cht$month == "November"
] <-
  calendar_days_2011_2019Cht$DayObserved[
    calendar_days_2011_2019Cht$year == 2017 &
      calendar_days_2011_2019Cht$month == "November"
  ] -2


calendar_days_2011_2019Cht$propMidyear<-(calendar_days_2011_2019Cht$DayObserved/calendar_days_2011_2019Cht$days_in_year)
calendar_days_2011_2019Cht<-calendar_days_2011_2019Cht%>%dplyr::select(year,month,propMidyear)

calendar_days_2011_2019Cht <- calendar_days_2011_2019Cht %>%
  dplyr::rename(
    yoa = year,
    admission_month = month
  ) %>%
  dplyr::mutate(
    yoa = as.numeric(yoa),
    admission_month = match(admission_month, month.name)
  )

#April    August  December  February   January      July      June     March       May  November 
#October September 
##readmore midyear population data and preparation
temp14frematunder1ALLn<-read_dta("data/temp14frematunder1ALL.dta")
temp14frematunder1ALLn$PersonMonthunder1_All<-((temp14frematunder1ALLn$midyrpopunder1)/12)
temp14frematunder1ALLn$PersonMonth5to14_All<-((temp14frematunder1ALLn$midyrpop5to14)/12)
temp14frematunder1ALLn<-temp14frematunder1ALLn%>%rename(midyrpopunder1_All=midyrpopunder1,midyrpop5to14_All=midyrpop5to14)

AdmissionInc<-inner_join(AdmissionInc, temp14frematunder1ALLn, by =c("yoa"))
AdmissionInc<-inner_join(AdmissionInc, calendar_days_2011_2019Cht, by =c("yoa","admission_month"))


##incidentDenominator(Try this please.I think it can be the best)
AdmissionInc$personMonthunder1_trybestIncDenominator<-(AdmissionInc$midyrpopunder1_All)*(AdmissionInc$propMidyear)
AdmissionInc$personMonth5to14_trybestIncDenominator<-(AdmissionInc$midyrpop5to14_All)*(AdmissionInc$propMidyear)


##Lets prepare infant mortality count,Birth weight,and number of live births.
## import live births
Live_Births<- read_excel("data/KEMRIDATA_BEN/Mid_year_pop_2011_2024.xlsx",sheet = "liveBirthsCorrect")
names(Live_Births)
Live_Births<- Live_Births %>%
  dplyr::mutate(
    yoa = as.integer(str_extract(Birth_year_month, "^[0-9]{4}")),
    admission_month = as.integer(str_extract(Birth_year_month, "(?<=m)[0-9]+"))
  )

Live_Births<-Live_Births%>%dplyr::select(-Birth_year_month)%>%dplyr::rename(LiveBirths=Numbers)
Live_Births<-Live_Births[Live_Births$yoa<=2019,]
Live_Births<-Live_Births[!(Live_Births$yoa==2011 & Live_Births$admission_month<4),]

##check BenardNew births data
library(haven)
library(dplyr)
library(lubridate)
FINAL_DATA_08_Apr_2026<- read_dta("data/NEWFINAL_DATA_18_May_2026.dta")
names(FINAL_DATA_08_Apr_2026)

attr(FINAL_DATA_08_Apr_2026$new_Location_birth, "labels")
class(FINAL_DATA_08_Apr_2026$new_Location_birth)
##current live births

LiveBirths_currentdata <- FINAL_DATA_08_Apr_2026 %>%
  mutate(
    yoa = year(child_dob),
    admission_month = month(child_dob)
  ) %>%
  group_by(yoa, admission_month) %>%
  summarise(
    LiveBirths_currentdata = n(),
    .groups = "drop"
  ) %>%
  arrange(yoa, admission_month)

#View(LiveBirths_currentdata)


#check missing
#Gestational age (I have decided am not using this for imputation)
tab <- table(FINAL_DATA_08_Apr_2026$gestation, useNA = "ifany")
prop.table(tab)

#Maternal age (Using this in the imputation model)
tab <- table(FINAL_DATA_08_Apr_2026$maternal_age_yrs, useNA = "ifany")
prop.table(tab) 
summary(FINAL_DATA_08_Apr_2026$maternal_age_yrs)
 #birth order (Not using in imputation model)
tab <- table(FINAL_DATA_08_Apr_2026$birth_order, useNA = "ifany")
prop.table(tab)

#Maternal education (Not using in the imputation model 47% missing)
tab <- table(FINAL_DATA_08_Apr_2026$maternal_education, useNA = "ifany")
prop.table(tab)  


#Location of  birth (In hospital or outside the hospital;I will use this in the imputation model)
tab <- table(FINAL_DATA_08_Apr_2026$new_Location_birth, useNA = "ifany")
prop.table(tab)

FINAL_DATA_08_Apr_2026$new_Location_birth <- factor(
  FINAL_DATA_08_Apr_2026$new_Location_birth
)

is.factor(FINAL_DATA_08_Apr_2026$new_Location_birth)

levels(FINAL_DATA_08_Apr_2026$new_Location_birth)





#Child year birth (I will use this in imputation model)
tab <- table(FINAL_DATA_08_Apr_2026$child_year_birth, useNA = "ifany")
prop.table(tab)
#child_month birth (Using this in imputation model)
tab <- table(FINAL_DATA_08_Apr_2026$Child_month_birth, useNA = "ifany")
prop.table(tab)


FINAL_DATA_08_Apr_2026 <- FINAL_DATA_08_Apr_2026 %>%
  mutate(
    Child_month_birthFactor = factor(Child_month_birth)
  )

is.factor(FINAL_DATA_08_Apr_2026$Child_month_birthFactor)

table(FINAL_DATA_08_Apr_2026$Child_month_birthFactor, useNA = "ifany")


##Shortest distance to KCH
tab <- table(FINAL_DATA_08_Apr_2026$dist_kchKm, useNA = "ifany")
prop.table(tab)





  
   #Birth_weight (79% missing)
   tab <- table(FINAL_DATA_08_Apr_2026$birth_weight, useNA = "ifany")
   prop.table(tab)
   summary(FINAL_DATA_08_Apr_2026$birth_weight)
   
   nrow(FINAL_DATA_08_Apr_2026[is.na(FINAL_DATA_08_Apr_2026$birth_weight),])
    #prop missing
   print(round((nrow(FINAL_DATA_08_Apr_2026[is.na(FINAL_DATA_08_Apr_2026$birth_weight), ]) /
                  nrow(FINAL_DATA_08_Apr_2026)) * 100, 2))
   
   
   FINAL_DATA_08_Apr_2026_u<-FINAL_DATA_08_Apr_2026%>%dplyr::select(birth_weight,dist_kchKm,
                                                           child_year_birth,Child_month_birthFactor,Child_month_birth,child_dob,
                                                           new_Location_birth,maternal_age_yrs)
   
   saveRDS(FINAL_DATA_08_Apr_2026_u,"data/FINAL_DATA_08_Apr_2026_u.rds")

   
   
##Merge Live_Births and admissions data
AdmissionInc<-inner_join(AdmissionInc, Live_Births, by =c("yoa","admission_month"))
AdmissionInc<-inner_join(AdmissionInc, LiveBirths_currentdata, by =c("yoa","admission_month"))

##Deaths in under1 and 5 to 14
Deaths_2011_2024 <- read_excel("data/KEMRIDATA_BEN/Deaths_2011_2024.xlsx")
Deaths_2011_2024<- Deaths_2011_2024 %>%
  dplyr::mutate(
    Dob = as.Date(Dob),
    dod = as.Date(dod)
  )
Deaths_2011_2024<-Deaths_2011_2024%>%
  dplyr::mutate(
    yoa= year(dod),
    admission_month= month(dod)
  )
deaths<-Deaths_2011_2024%>%
  dplyr::filter(
    dod >= as.Date("2011-04-01") &
      dod <= as.Date("2019-12-31")
  )
deaths<-deaths%>%dplyr::arrange(dod)
view(deaths)

deaths<-deaths%>%dplyr::mutate(deathsunder1=1)
deaths <- deaths %>%
  dplyr::mutate(
    age_at_death1=(as.numeric(dod - Dob) /365.25),
    age_at_death2=floor(age_at_death1),
    age_months =floor(as.numeric(dod - Dob) / (365.25 / 12))
  )

#deaths <- deaths %>%
  #mutate(
    #age_months = as.numeric(dod - Dob) / (365.25 / 12)
  #)

all.equal(deaths$age_at_death2, deaths$age_at_death)
deaths<-deaths%>%dplyr::select(-age_at_death1,-age_at_death,-age_months)

deaths <- deaths %>%
  dplyr::mutate(
    age_category = case_when(
      age_at_death2 < 1 ~ "<1 year",
      age_at_death2 >= 5 & age_at_death2 < 15 ~ "5–14 years",
      TRUE ~ NA_character_
    )
  )
deaths<-deaths[deaths$age_at_death2>=0,]

deathsgroup<-deaths%>%group_by(yoa, admission_month,age_category) %>%
  summarise(
    deathsunder1 = sum(deathsunder1, na.rm = TRUE),
    .groups = "drop"
  )

deaths_clean <- deathsgroup %>%
  dplyr::filter(!is.na(age_category))

deaths_wide <- deaths_clean %>%
  pivot_wider(
    id_cols = c(yoa, admission_month),
    names_from = age_category,
    values_from = deathsunder1,
    values_fill = 0
  )

view(deaths_wide)

deaths_wide<-deaths_wide%>%dplyr::rename(deaths5to14=`5–14 years`,deathsunder1=`<1 year`)

##Merge admission and deaths
AdmissionInc<-inner_join(AdmissionInc, deaths_wide, by =c("yoa","admission_month"))
AdmissionInc<-AdmissionInc%>%dplyr::rename(midyrpop_under1=midyrpopunder1,person_Month_under1=person_Month)
view(AdmissionInc)
##Lets generate the incidences
names(AdmissionInc)

AdmissionInc$personMonthunder1_trybestIncDenominatorPLus1<-(AdmissionInc$personMonthunder1_trybestIncDenominator)+1
AdmissionInc$personMonth5to14_trybestIncDenominatorPLus1<-(AdmissionInc$personMonth5to14_trybestIncDenominator)+1


AdmissionInc$Malaria_Incidence<-((AdmissionInc$Malaria_count/AdmissionInc$person_Month_under1)*100000)
AdmissionInc$Malaria_Incidenceb<-((AdmissionInc$Malaria_count/AdmissionInc$personMonthunder1_trybestIncDenominator)*100000)
AdmissionInc$Malaria_count<-NULL

AdmissionInc$Measles_Incidence<-((AdmissionInc$Measles_count/AdmissionInc$person_Month_under1)*100000)
AdmissionInc$Measles_Incidenceb<-((AdmissionInc$Measles_count/AdmissionInc$personMonthunder1_trybestIncDenominator)*100000)
AdmissionInc$Measles_count<-NULL

AdmissionInc$birth_defects_Incidence<-((AdmissionInc$birth_defects_count/AdmissionInc$person_Month_under1)*100000)
AdmissionInc$birth_defects_Incidenceb<-((AdmissionInc$birth_defects_count/AdmissionInc$personMonthunder1_trybestIncDenominator)*100000)
AdmissionInc$birth_defects_count<-NULL

AdmissionInc$malnutrition_Incidence<-((AdmissionInc$malnutrition_count/AdmissionInc$person_Month_under1)*100000)
AdmissionInc$malnutrition_Incidenceb<-((AdmissionInc$malnutrition_count/AdmissionInc$personMonthunder1_trybestIncDenominator)*100000)
AdmissionInc$malnutrition_count<-NULL

AdmissionInc$HIV_Incidence<-((AdmissionInc$HIV_count/AdmissionInc$person_Month_under1)*100000)
AdmissionInc$HIV_Incidenceb<-((AdmissionInc$HIV_count/AdmissionInc$personMonthunder1_trybestIncDenominator)*100000)
AdmissionInc$HIV_count<-NULL

AdmissionInc$neonatal_sepsis_Incidence<-((AdmissionInc$neonatal_sepsis_count/AdmissionInc$person_Month_under1)*100000)
AdmissionInc$neonatal_sepsis_Incidenceb<-((AdmissionInc$neonatal_sepsis_count/AdmissionInc$personMonthunder1_trybestIncDenominator)*100000)
AdmissionInc$neonatal_sepsis_count<-NULL

AdmissionInc$birth_Asphyxia_Incidence<-((AdmissionInc$birth_Asphyxia_count/AdmissionInc$person_Month_under1)*100000)
AdmissionInc$birth_Asphyxia_Incidenceb<-((AdmissionInc$birth_Asphyxia_count/AdmissionInc$personMonthunder1_trybestIncDenominator)*100000)
AdmissionInc$birth_Asphyxia_count<-NULL

saveRDS(AdmissionInc,"data/AdmissionInc.rds")


##generate time variable
AdmissionInc <- AdmissionInc %>%
  dplyr::arrange(yoa, admission_month) %>%
  dplyr::mutate(time = row_number())
##Policy variable
AdmissionInc <- AdmissionInc %>%
  dplyr::mutate(
    Policy_Indicator= if_else(
      yoa > 2013 | (yoa == 2013 & admission_month >= 7),
      1, 0
    )
  )
view(AdmissionInc)
AdmissionInc<-AdmissionInc%>%dplyr::mutate(Policy_Indicator=ifelse((yoa == 2013 & admission_month==6),NA,Policy_Indicator))

##March 2017 time 72
AdmissionInc <- AdmissionInc %>%dplyr::mutate(Policy_Interraction=((time-72)*Policy_Indicator))
##factor month
#delta <- 0.0001
AdmissionInc <- AdmissionInc %>%
  dplyr::arrange(time) %>%
  dplyr::mutate(
    month_dummy= factor(admission_month, levels = 1:12, labels = month.abb),
    #t_c   = t - mean(t, na.rm = TRUE),
    #off   = log(control + delta)     # <-- offset with small positive delta
  )



##checking if control is nice(My control real seem to be working.See unaffected by policy and meeting common trend assumption)
## Pre-policy data
df_pre <- AdmissionInc %>%
  dplyr::filter(Policy_Indicator == 0)

## Fit NB model
m_nb <- MASS::glm.nb(
  deathsunder1 ~ time + month_dummy + offset(log(deaths5to14)),
  data = df_pre
)

## Newey–West lag
n <- nrow(df_pre)
L <- floor(n^(1/4))

## Newey–West variance–covariance matrix
V_hac <- sandwich::NeweyWest(
  m_nb,
  lag = L,
  prewhite = FALSE,
  adjust = TRUE
)

## ---- Bottomley-style Wald test ----
beta_time <- coef(m_nb)["time"]
se_time   <- sqrt(diag(V_hac))["time"]

z_time <- beta_time / se_time
p_time <- 2 * (1 - pnorm(abs(z_time)))

## Results
beta_time
se_time
z_time
p_time

##Is my control affected by the policy;
m_ctrl <- MASS::glm.nb(
  deaths5to14 ~ time+ Policy_Indicator+ month_dummy,
  data =AdmissionInc
)

