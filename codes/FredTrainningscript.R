setwd("C:/Users/THINKPAD/Desktop/PleliminaryAnalysis")
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
