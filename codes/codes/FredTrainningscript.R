setwd("C:/Users/THINKPAD/Desktop/PleliminaryAnalysis")
##packages
require(haven)
require(dplyr)
require(tidyr)
admin <- read_dta("Data/temp/admin.dta")##import data
##generate resident status
admin$resident<-NULL
admin$resident<-NA
admin$resident[!is.na(admin$fk_res)]<-1
admin$resident[is.na(admin$fk_res)]<-0
admin$resident<-factor(admin$resident,levels=c(1,0),labels =c("dss","non-dss") )

##view resident
View(admin%>%select(resident,fk_res,everything()))
##Table
table(admin$resident)
##generating year of study and month amd week
admin<-admin%>%mutate(age_days=(doa-dob))
View(admin%>%select(age_days,doa,dob,everything()))
class(admin$age_days)
admin$age_days<-as.numeric(admin$age_days)
#age in months
admin$agemonths<-round(admin$age_days/30.4375)
admin$agemonths<-as.numeric(admin$agemonths)
View(admin%>%select(age_days,agemonths,doa,dob,everything()))
#3Age in years
admin$ageyears<-round(admin$age_days/365.25)
admin$ageyears<-as.numeric(admin$ageyears)
##Age category
admin$agecat<-ifelse(admin$ageyears<5,0,ifelse(admin$ageyears>=5,1,NA))
admin$agecat<-factor(admin$agecat,levels =c(0,1),labels=c("<5yrs",">=5yrs") )
View(admin%>%select(age_days,agemonths,ageyears,agecat,doa,dob,everything()))
##*drop if age_days<0 // 5 records droppped
admin<-admin%>%filter(ageyears>=0 & ageyears<=14)
min(admin$ageyears)
max(admin$ageyears)
##YEAR,MONTH AND WEEK
admin$yearadmission<-substr(as.POSIXct(admin$doa),1,4)
admin$monthadmision<-months(admin$doa)
admin$weekadmission<-weekdays(admin$doa)
View(admin%>%select(age_days,agemonths,ageyears,agecat,yearadmission,monthadmision,weekadmission,doa,dob,everything()))
#Duplicates removal
nrow(admin)
admin<-admin%>%distinct()
admin<-distinct(admin)
##For duplicates based on columns remember to pt the option .Keep_all=TRUE so that R doesnt drop other variables
admin<-admin%>%distinct(fk_person,doa,.keep_all =TRUE)
View(admin)
##Library lubridate or generating year months day time
require(lubridate)
admin$yearlub<-year(admin$doa)
admin$monlub<-month(admin$doa)
admin$daylub<-day(admin$doa)
##time 
admin$toa<-hms(admin$time_admn)
admin$hour<-hour(admin$toa)
admin$minutes<-minute(admin$toa)
admin$seconds<-second(admin$toa)

subset(admin,select=c(time_admn,toa,hour,minutes,seconds))






