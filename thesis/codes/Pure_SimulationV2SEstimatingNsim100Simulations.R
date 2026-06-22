
##Package area.Load packages
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
library(tscount)
#setting seed
set.seed(123)
#parameters
x1drop_schock<--0.9 
driftX1<-0.000 #0.015
seasonalX<-c(log(1.04),log(0.98),log(1.04),log(1.01))
lagX2<-4
lagtq<-2
lagpstT0X1<-15
#sim_listUncommon<-list()
vcoeffsY<-c(log(700),log(0.995),log(0.96),(log(0.84)/(x1drop_schock)),350)#vary B2 (-0.04,-0.36,-0.51) #B0,B1(time),B2(policy),B3(confounder),thetaY repectilvely
vcoeffsZ<-c(log(900),log(0.995),(log(0.84)/(x1drop_schock)),350) # bo,b1(time),b2(confounder),thetaZ.Lets tKE 90% trend,coef
EstimandMain<- log(0.70)
fredSiProject<-function(nsim,seasonalX,vcoeffsY,vcoeffsZ,rhofd, sigma2fd,n,u_s=0){ ##Let var be 3 times the mean. mean to variance ration for choosing overdispersion.Chose mean to var ratio of 2;mean=B0
sim_listcommon<-list()
  for (i in 1:nsim){
    ##read coefficients
    ##Y
    B0<-vcoeffsY[1]
    B1<-vcoeffsY[2]
    B2<-vcoeffsY[3]
    B3<-vcoeffsY[4]
    thetaY<-vcoeffsY[5]
    ##Cofficient Z
    b0<-vcoeffsZ[1]
    b1<-vcoeffsZ[2]
    b2<-vcoeffsZ[3]
    thetaZ<-vcoeffsZ[4]
    
    
    
    T<-n #100 for now
    t  <- 1:T
    #t0 <- 25
    if (T%%2==0){
      t0<-((T/2)+1)
      
    }else{
      t0<- floor((T/2)+0.5) 
    }
    
    print(t0)
    
    P  <- as.integer(t >= t0)     # policy step at t=25
    #N  <- rep(1000, T)            # exposure (can vary if you want)
    
    ##letting confounder have level drop 3 step early around policy
    tq<-t0-lagtq
    print(tq)
    ## ----- CONFOUNDERS (plausible, not collinear with policy) -----
    # X1: near-coincident shock at t>=22 + small drift + noise
    #seasonalX<-c(log(1.04),log(0.98),log(1.04),log(1.01)) #Just picked berneal et al to inform this seasonality component
    Season<-seasonalX[1]*sin((2*pi*t*1)/12)+seasonalX[2]*cos((2*pi*t*1)/12)+seasonalX[3]*sin((2*pi*t*2)/12)+seasonalX[4]*cos((2*pi*t*2)/12)
    X1 <- x1drop_schock* as.integer(t >=tq & (t<=(t0+lagpstT0X1))) + driftX1* (t)+Season + rnorm(T, 0, 0.20)
    X1 <- X1 - mean(X1[t < t0])
    #X1<- as.integer(t >= (t0 - 3))
    # Seasonality (optional but realistic)
    #S1 <- sin(2*pi*t/12)
    #C1 <- cos(2*pi*t/12)
    
    ## Quick checks for (non-)collinearity with policy
    cor_P_X1 <- cor(P, X1)   # should be far from 1
    print(round(cor_P_X1, 3))
    
    ## ----- TRUE COEFFICIENTS -----internd to vary B2(-0.36 expected,-0.51 maximum,-0.04 minimal)
    # Outcome Y (your estimand is B2)
    #B0 <-3.99  #Bo infomred by mean pneumonia counts 2002
    #B1 <- -0.005
    #B2 <- -0.36     # true policy effect on log-mean for Y (we want to recover this)
    #g1 <- ln(1.75)/(-0.9)
    #B3<--0.623    # effect of confounder X1 on Y(Informed by strike effect Ongayo)
    #g2a <- 0.25     # seasonality on Y
    #g2b <- -0.15
    #thetaY <- 20    # NB2 dispersion (larger -> less overdispersion)
    
    # Control Z: similar trend and confounders, but NO policy effect(I want to above the outcome slightly)
    #b0 <-4.15  #Bo infomred by mean pneumonia counts 2002
    #b1 <- -0.005  #-0.004
    #b2<--0.623  # -0.073
    #g2az <- 0.20
    #g2bz <- -0.10
    #thetaZ <- 25
    
    
    ##Simulate auto-corellated errors
    ## --- AR(1) error for the log-mean ---
    rho    <-rhofd        #0.4 # choose from {0, 0.2, 0.4, 0.6, 0.8}/I took the mean Turner
    sigma2 <-sigma2fd         #0.1 # variance of white-noise w_t (per Turner et al.)
    
    u <- numeric(T)
    sd_stat <- sqrt(sigma2 / (1 - rho^2))  # stationary SD of AR(1) error
    u[1] <- rnorm(1, 0, sd_stat)
    for (tt in 2:T) u[tt] <- rho * u[tt-1] + rnorm(1, 0, sqrt(sigma2))
    
    ## Optional (recommended for counts): mean-preserve on log scale so E[exp(u)] ~ 1
    u <- u - 0.5 * (sigma2 / (1 - rho^2)) #ensure
    
    if (u_s==1){
      ## ----- GENERATE COUNTS (NB2) -----
      # Linear predictors (log-means)
      ##common trend data
      etaY <- B0 + B1 * t + B2 * P + B3 * X1+u
      muY  <- exp(etaY)
      Y    <- rnbinom(T, size = thetaY, mu = muY)
      
      etaZ <- b0 + b1 * t+ b2* X1+u
      muZ  <- exp(etaZ)
      Z    <- rnbinom(T, size = thetaZ, mu = muZ)
      
      ##unshared Z confounder
      etaZ <- b0 + b1 * t+u
      muZ  <- exp(etaZ)
      Z_unXrm<- rnbinom(T, size = thetaZ, mu = muZ)
      
      dat<- data.frame(t, P,X1,Y, Z, Z_unXrm,t0=t0,tq=tq)
      dat$j<-i
      sim_listcommon[[i]]<-dat
      ##unshared trend data
      #etaY <- B0 + B1 * t + B2 * P + g1 * X1+u
      #muY  <- exp(etaY)
      #Y    <- rnbinom(T, size = thetaY, mu = muY)
      
     
      
    }else{
      ## ----- GENERATE COUNTS (NB2) -----
      # Linear predictors (log-means)
      ##common trend data
      etaY <- B0 + B1 * t + B2 * P + B3 * X1
      muY  <- exp(etaY)
      Y    <- rnbinom(T, size = thetaY, mu = muY)
      
      etaZ <- b0 + b1 * t+ b2* X1
      muZ  <- exp(etaZ)
      Z    <- rnbinom(T, size = thetaZ, mu = muZ)
      
      ##unshared Z confounder
      etaZ <- b0 + b1 * t
      muZ  <- exp(etaZ)
      Z_unXrm<- rnbinom(T, size = thetaZ, mu = muZ)
      
      dat<- data.frame(t, P,X1,Y,Z,Z_unXrm,t0=t0,tq=tq)
      dat$j<-i
      sim_listcommon[[i]]<-dat
    }
    
  }
  ## return BOTH
  #list(common = sim_listcommon, uncommon = sim_listUncommon)
  list(common = sim_listcommon)
  }
  
##running this function to generate shared and unshared confounder# run u_s=0 for not adding AR(1)
out <-fredSiProject(100,seasonalX,vcoeffsY,vcoeffsZ,0.8,0.1,12,u_s =1)
names(out) 
#d1_common   <- out$common[[1]]
#d1_uncommon <- out$uncommon[[1]]
##Bind all shared confounder simulation
require(dplyr)
common_all <- do.call(rbind, out$common)
common_all<-common_all%>%arrange(j,t0)
t0<-common_all$t0[1]
tq<-common_all$tq[1]
print(t0)
print(tq)
##drop these unnecessary columns
common_all$t0<-NULL
common_all$tq<-NULL

common_all$tce<-common_all$t-t0
common_all$Policy_Time<-((common_all$tce)*(common_all$P))
View(common_all)
##END OF SIMULATION




##ANALYSIS
##plot simulation
## Pre-policy trends overtime.
#Shared confounder common trend loess using first simulated data;
dat<-common_all[common_all$j==1,]
datpre <- subset(dat, P == 0|P == 1)

datpre <- subset(dat, P == 0|P==1)

##check corelation
cor(dat[, c("t","tce","P","X1","Y","Z")], use = "pairwise.complete.obs", method = "pearson")


fitY <- loess(Y ~ t, data = datpre,span=0.6, degree = 1, family = "symmetric")
fitZ <- loess(Z ~ t, data = datpre, span=0.6,degree = 1, family = "symmetric")

plot(datpre$t, datpre$Y, pch = 16, col = "steelblue",
     xlab = "Time (pre-policy)", ylab = "Count",
     main = "Estimated trend(A)")
points(datpre$t, datpre$Z, pch = 16, col = "tomato")
lines(datpre$t, predict(fitY), lwd = 3, col = "steelblue4")
lines(datpre$t, predict(fitZ), lwd = 3, col = "tomato4")
abline(v=t0,    lty =8, col = "black") 
grid()

legend("topleft",
       legend = c("Y LOESS", "Z LOESS"),
       col    = c("steelblue4", "tomato4"),
       lty    = c(1, 1),
       lwd    = c(3, 3),
       bty = "n", inset = 0.01)
##common trend model
##testing common trend first
datpre<-subset(dat,P==0)
m_off1<- glm.nb(Y ~t + offset(log(Z)),
                data =datpre)
summary(m_off1)

## Install + load HAC tools
if (!requireNamespace("sandwich", quietly = TRUE)) install.packages("sandwich")
if (!requireNamespace("lmtest",    quietly = TRUE)) install.packages("lmtest")
library(sandwich)
library(lmtest)

L <-floor(nrow(dat)^(1/4))  # Wooldridge
print(L)
# Compute Newey-West variance-covariance matrix (lag = 3)
vcov_nw <- NeweyWest(m_off1, lag =L, prewhite = FALSE, adjust = TRUE)

# Print coefficient table with Newey-West standard errors
coeftest(m_off1, vcov = vcov_nw)






## ----- FIT MODELS -----
# We'll use MASS::glm.nb for NB2
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)


##generate storing variables
common_all$dataset<-NA

#confimr observations
common_all$Allobs<-NA
common_all$AllobsPre<-NA
#extract p_valuefor common tren test
common_all$CommonAsssumption_Pvalue<-NA

##CITS
common_all$P_estimateCITS<-NA
common_all$Pstd_estimateCITS<-NA

##Traditional regression model adjusted
common_all$P_estimateTrd<-NA
common_all$Pstd_estimateTrd<-NA
gh<-max(common_all$j)
print(gh)
kl<-1

##cheking in how many simulations do I meet common trend assumption Assumption?
countPlessAlphaCITS=0
##checking how many of My Z are zero
sumzeroZ=0
sumnonzeroZ=0
common_all<-common_all%>%arrange(j,t)
for (prof in 1:gh) {
  common_all<-common_all%>%arrange(j,t)
  common_all$dataset[kl]<-prof
  nk<-nrow(common_all[common_all$j==prof,])
  nkpre<-nrow(common_all[((common_all$j==prof) & (common_all$t<t0)),])
  common_all$Allobs[kl]<-nk
  common_all$AllobsPre[kl]<-nkpre
  nkkdt<-common_all[common_all$j==prof,]
  klf<-min(nkkdt$Z)
  if (klf<=0){
    print(paste("ZERO is in dataset",prof,"for Z",sep =""))
    sumzeroZ<-sumzeroZ+1 
  }else{
    sumnonzeroZ<-sumnonzeroZ+1
  }
  ##testing common trend first
  #datpre<-subset(dat,P==0)
  prenkt<-common_all[((common_all$j==prof) & (common_all$t<t0)),]
  m_off1<- glm.nb(Y ~t + offset(log(Z)),
                  data =prenkt)
  
  # Compute Newey-West variance-covariance matrix (lag = 3)
  vcov_nw <- NeweyWest(m_off1, lag =L, prewhite = FALSE, adjust = TRUE)
  
  # Print coefficient table with Newey-West standard errors
  ct11<-coeftest(m_off1, vcov = vcov_nw)
  
  koblo<- ct11["t", grep("^Pr\\(", colnames(ct11), value = TRUE)]
  
  
  ##extracting Evidence of meeting common trend assumption
  common_all$CommonAsssumption_Pvalue[kl]<-koblo
  if (koblo<0.05){
    countPlessAlphaCITS=countPlessAlphaCITS+1 
  }else{
    countPlessAlphaCITS<-countPlessAlphaCITS+0
  }
  
  
  
  #print(nk)
  ##fitting common trend model with control as offset
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
  m_off <- glm.nb(Y ~P+offset(log(Z)),
                  data =nkkdt)
  # Compute Newey-West variance-covariance matrix (lag = 3)
  vcov_nw <- NeweyWest(m_off, lag =L, prewhite = FALSE, adjust = TRUE)
  
  # Print coefficient table with Newey-West standard errors
  ct<-coeftest(m_off, vcov = vcov_nw)
  
  ##fill coefficeints
  common_all$P_estimateCITS[kl]<-unname(ct["P","Estimate"])
  common_all$Pstd_estimateCITS[kl]<-unname(ct["P","Std. Error"])
  
  
  
  
  ##traditional regression 
  ##with uncentered time
  m_adj <- glm.nb(Y ~ t + P + X1, data =nkkdt)
  # Compute Newey-West variance-covariance matrix (lag = 3)
  vcov_nw <- NeweyWest(m_adj, lag =L, prewhite = FALSE, adjust = TRUE)
  
  # Print coefficient table with Newey-West standard errors
  trr<-coeftest(m_adj, vcov = vcov_nw)
  #fill coefficients
  common_all$P_estimateTrd[kl]<-unname(trr["P","Estimate"])
  common_all$Pstd_estimateTrd[kl]<-unname(trr["P","Std. Error"])
  
  kl<-kl+1
}
##print number of zeroes in Z if any
print(sumzeroZ)
print(sumnonzeroZ)
print(countPlessAlphaCITS)
##mean of CITS P and standard error
##mean coef
CITS_Pmean<-mean(common_all$P_estimateCITS,na.rm =TRUE)
print(CITS_Pmean)
##mean std
CITS_Pstdmean<-mean(common_all$Pstd_estimateCITS,na.rm =TRUE)
print(CITS_Pstdmean)


##mean traditional
#mean coef
Trd_Pmean<-mean(common_all$P_estimateTrd,na.rm=TRUE)
print(Trd_Pmean)
Trd_Pstdmean<-mean(common_all$Pstd_estimateTrd,na.rm=TRUE)
print(Trd_Pstdmean)


##Lets Add columns for estimating nsmim for bias and MSE(ONLY FOR CALCULATING Nsim) #Used 100 sim and 150 obs(T=150) for this nsim estimation
## seed 123(Tre)
##Nsim Bia (Trd=2435;CITS=156)
common_all$ThetaHatMINUS_ThetaTrueCITS<-(common_all$P_estimateCITS-CITS_Pmean)
common_all$ThetaHatMINUS_ThetaTrueTrd<-(common_all$P_estimateTrd-Trd_Pmean)
BiasNsimTrdCITS<-((sum((common_all$ThetaHatMINUS_ThetaTrueCITS)^2, na.rm = TRUE)/(gh-1))/(0.005^2))
print(ceiling(BiasNsimTrdCITS)) ##both give nsim required given CITS
print(ceiling(var(common_all$P_estimateCITS,na.rm =TRUE)/(0.005)^2)) ##both give nsim required given CITS
#((sum((common_all$ThetaHatMINUS_ThetaTrueTrd)^2, na.rm = TRUE)/(gh-1))/(0.005^2))
BiasNsimTrd<-((sum((common_all$ThetaHatMINUS_ThetaTrueTrd)^2, na.rm = TRUE)/(gh-1))/(0.005^2))
print(ceiling(BiasNsimTrd)) ##both give nsim required given Trd
print(ceiling(var(common_all$ThetaHatMINUS_ThetaTrueTrd,na.rm =TRUE)/(0.005)^2)) ##both give nsim required given Trd

##MSE(Nsim)
##CITS
EstimandMain<-log(0.96)
xmse<-(common_all$P_estimateCITS-EstimandMain)^2
ymse<-(var(xmse,na.rm = TRUE)/(0.005)^2)
print(ceiling(ymse))
MSECITS<-sum((common_all$P_estimateCITS-EstimandMain)^2,na.rm=TRUE)/(gh)
nsimMSCECITS<-(sum(((common_all$P_estimateCITS-(EstimandMain))^2-MSECITS)^2,na.rm =TRUE)/(gh-1))/(0.005^2)
print(nsimMSCECITS)  ##Nsim required CITS
##Traditional regression
EstimandMain<-log(0.96)
xmse<-(common_all$P_estimateTrd-EstimandMain)^2
ymse<-(var(xmse,na.rm = TRUE)/(0.005)^2)
print(ceiling(ymse))
MSETrd<-sum((common_all$P_estimateTrd-EstimandMain)^2,na.rm=TRUE)/(gh)
nsimMSCETrd<-(sum(((common_all$P_estimateTrd-(EstimandMain))^2-MSETrd)^2,na.rm =TRUE)/(gh-1))/(0.005^2)
print(nsimMSCETrd)  ##Nsim required TRd

##Nsim for EmpericalSE
##trd
varem<-var(common_all$P_estimateTrd,na.rm =TRUE)
morva<-((varem)/(2*(0.005)^2))
nsimtrdEMPSE<-(1+morva)
ceiling(print(nsimtrdEMPSE)) #sim based on trd
#cits
varem<-var(common_all$P_estimateCITS,na.rm =TRUE)
morva<-((varem)/(2*(0.005)^2))
nsimtrdEMPSE<-(1+morva)
ceiling(print(nsimtrdEMPSE)) #sim based on trd

##Nsim for Modelbased MSE
## ---------- Nsim for Model-based SE (Trd and CITS) ----------

target_mcse_ModSE <- 0.005   # your chosen MCSE target on log scale

## --- Traditional regression (Trd) ---
# model-based SEs from each simulation:
se_trd  <- common_all$Pstd_estimateTrd

# average model-based SE (ModSE_hat):
ModSE_trd <- sqrt(mean(se_trd^2, na.rm = TRUE))

# variance of the variance estimates Var_hat(theta) ≈ Var(se^2):
vhat_trd      <- se_trd^2
var_varhat_trd <- var(vhat_trd, na.rm = TRUE)

# required n_sim for model-based SE:
nsim_ModSE_trd <- var_varhat_trd /
  (4 * ModSE_trd^2 * target_mcse_ModSE^2)

ceiling(nsim_ModSE_trd)  # N_sim based on Trd


## --- Controlled Interrupted Time Series (CITS) ---
se_cits  <- common_all$Pstd_estimateCITS

ModSE_cits <- sqrt(mean(se_cits^2, na.rm = TRUE))

vhat_cits       <- se_cits^2
var_varhat_cits <- var(vhat_cits, na.rm = TRUE)

nsim_ModSE_cits <- var_varhat_cits /
  (4 * ModSE_cits^2 * target_mcse_ModSE^2)

ceiling(nsim_ModSE_cits)  # N_sim based on CITS




# Load library
library(ggplot2)

# Create data
gbv_data <- data.frame(
  Category = c("Reported GBV", "No GBV Reported", "No Response"),
  Percentage = c(39, 50, 11)
)

# Bar plot
ggplot(gbv_data, aes(x = Category, y = Percentage)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(Percentage, "%")),
            vjust = -0.5, size = 5) +
  ylim(0, 60) +
  labs(
    title = "Household Reports of Gender-Based Violence (GBV)",
    x = "",
    y = "Percentage of Households"
  ) +
  theme_minimal(base_size = 14)




