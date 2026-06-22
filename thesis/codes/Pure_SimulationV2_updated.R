
##Package area.Load packages
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

#if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)

#if (!requireNamespace("withr", quietly = TRUE)) install.packages("withr")
library(withr)

#if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
#library(tscount)
#parameters
x1drop_schock<--0.9 
driftX1<-0.000 #0.015
seasonalX<-c(log(1.04),log(0.98),log(1.04),log(1.01))
#Season<-
#harmonic(month, 2, 12)1  0.0370967  0.0100241    3.701 0.000520 ***
  #harmonic(month, 2, 12)2 -0.0182122  0.0096435   -1.889 0.064537 .  
#harmonic(month, 2, 12)3  0.0381608  0.0096848    3.940 0.000244 ***
  #harmonic(month, 2, 12)4  0.0147634  0.0097221    1.519 0.134935 

#sim_listUncommon<-list() choose bernal paper on Its to inform Bo,for Y 600,for z 700 #(log(1.75)/(-0.9))
vcoeffsY<-c(log(700),log(0.995),log(0.96),log(0.70),log(0.60),(log(0.84)/(x1drop_schock)),log(1),350)#vary B2 (-0.04,-0.36,-0.51) #B0,B1(time),B2(policy),B3(confounder),thetaY repectilvely
vcoeffsZ<-c(log(900),log(0.995),(log(0.84)/(x1drop_schock)),log(1.05),log(1),350) # bo,b1(time),b2(confounder),thetaZ.Lets tKE 90% trend,coef
#PolicyEffectOnZ<-c(log(0.985),log(0.97),log(0.95)) previous but let me make it 5% 10% 15%(also oposite for 15%)  too as below
PolicyEffectOnZ<-c(log(0.95),log(0.90),log(0.85),log(1.15))
#policy on X
#PolicyEffectOnX<-c(0.015,0.03,0.05)
PolicyEffectOnX<-c(log(1.05),log(1.10),log(1.15)) #let me try 5% 10% 15% increase on X
##X parameters
paraneterX<-c(log(2000),log(1.001),log(0.98))
lagX2<-4
lagtq<-2
lagpstT0X1<-15
##differenZslope
DifferentZSlope<-c(-0.0035,-0.0020,-0.0005,0.001)
## -0.0001 tune qudratic term its a problem
EstimandMain<- -0.36
#rhofd<-c(0,0.2,0.4,0.6,0.8)
rhofd<-c(0.4,0.2,0,0.6,0.8)
#{0, 0.2, 0.4, 0.6, 0.8}
#set seed
#RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
#set.seed(123)
n<-c(12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,80,100,150)
fredSiProject<-function(nsim,seasonalX,vcoeffsY,vcoeffsZ,rhofd,Rhosingle=0.4,SingleRho=FALSE, sigma2fd,nseries,u_s=0,genX2=TRUE,zslopPercent=10,PincZsl=10,averaging_n=150,shared_conf= 0.85,rho_z_const= 0.4,prop_MedX2=0.2,tensX2=1000,meanpreX2=2000,seed=123){ ##Let var be 3 times the mean. mean to variance ration for choosing overdispersion.Chose mean to var ratio of 2;mean=B0
  res <- list()
  pos0p4<-which(rhofd == 0.4)
  alpha<-shared_conf
  rho_z_const<-rho_z_const
  ##function for generating AR terms
  gen_y_twoz_ar1 <- function(T,
                             rho_y,
                             rho_z_var,
                             rho_z_const = 0.4,
                             sigma2,
                             alpha){
    
    # clamp alpha to valid range
    alpha <- max(-1, min(1, alpha))
    
    # 1) Innovation covariance matrix (same for both scenarios)
    Sigma_e <- sigma2 * matrix(c(1, alpha,
                                 alpha, 1),
                               nrow = 2, byrow = TRUE)
    
    # 2) Draw bivariate normal innovations: each marginal N(0, sigma2),
    #    Corr(e_y, e_z)=alpha
    E <- mvrnorm(n = T, mu = c(0, 0), Sigma = Sigma_e)
    e_y <- E[, 1]
    e_z <- E[, 2]   # this same e_z will drive both Z series
    
    # 3) Generate AR(1) processes
    u_y <- numeric(T)
    u_z_var <- numeric(T)
    u_z_const <- numeric(T)
    
    # Start from stationarity (marginal)
    u_y[1]       <- rnorm(1, 0, sqrt(sigma2 / (1 - rho_y^2)))
    u_z_var[1]   <- rnorm(1, 0, sqrt(sigma2 / (1 - rho_z_var^2)))
    u_z_const[1] <- rnorm(1, 0, sqrt(sigma2 / (1 - rho_z_const^2)))
    
    for(t in 2:T){
      u_y[t]       <- rho_y      * u_y[t-1]       + e_y[t]
      u_z_var[t]   <- rho_z_var  * u_z_var[t-1]   + e_z[t]
      u_z_const[t] <- rho_z_const* u_z_const[t-1] + e_z[t]
    }
    
    # ---- mean-preserving correction on log scale ----
    var_y       <- sigma2 / (1 - rho_y^2)
    var_z_var   <- sigma2 / (1 - rho_z_var^2)
    var_z_const <- sigma2 / (1 - rho_z_const^2)
    
    u_y       <- u_y       - 0.5 * var_y
    u_z_var   <- u_z_var   - 0.5 * var_z_var
    u_z_const <- u_z_const - 0.5 * var_z_const  
    
    list(u_yt =u_y,
         u_ct= u_z_var,
         u_ctconstant= u_z_const,
         e_y = e_y,
         e_z = e_z,
         params = list(rho_y=rho_y, rho_z_var=rho_z_var, rho_z_const=rho_z_const,
                       sigma2=sigma2, alpha=alpha))
  }
  
  #setting seed
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
#sim_listcommonBO_Expected<-list()
#sim_listcommonBO_MinEf<-list()
#sim_listcommonBO_MaxEf<-list()
sim_listcommonBO_MinEf <- list()
sim_listcommonBO_Expected <- list()
sim_listcommonBO_MaxEf <- list()
#sim_listcommon<-list()
u_ctconstant_store <- list()
for (rof in rhofd) {
  if (SingleRho){
    rho<-Rhosingle
    rof<-Rhosingle
  }else{
    rho<-rof
  }
  for(jojog in nseries){
    n<-jojog
    ##resetting seed to 0.4 rho jojog
    #if(rof==0.4 & n==averaging_n){
      #set.seed(seed)
    #}
    
    cat("corelation", rof,"of dataset",jojog,"(","18 diferent series lengths",")","\n")
    #cat("Length of series is", jojog, "\n")
    
    #rho<-rof
    oyaoya<-1
    for (jek in 3:5) {
      sim_listcommon<-list()
      ##read coefficients
      ##Y
      B0<-vcoeffsY[1]
      B1<-vcoeffsY[2]
      B2<-vcoeffsY[jek]
      B3<-vcoeffsY[6]
      B4<-vcoeffsY[7] #curve
      thetaY<-vcoeffsY[8]
      
      
      for (i in 1:nsim){
        #print(paste("simulating dataset",i,"Level change",B2,sep=""))
        
        ##Cofficient Z
        b0<-vcoeffsZ[1]
        b1<-vcoeffsZ[2]
        b2<-vcoeffsZ[3]
        b3<-vcoeffsZ[4]
        b4<-vcoeffsZ[5] #carvature
        thetaZ<-vcoeffsZ[6]
        #vcoeffsZ<-c(log(900),log(0.995),(log(0.84)/(x1drop_schock)),log(1.05),log(1),350) # bo,b1(time),b2(confounder),thetaZ.Lets tKE 90% trend,coef
        #Spillover:
        #PolicyEffectOnZ<-c(-0.015,-0.03,-0.05)
        #PolicyEffectOnZ<-c(log(0.95),log(0.90),log(0.85))
        #PolicyEffectOnZ<-c(log(0.95),log(0.90),log(0.85),log(1.15))
        b3_1Pfivepct<-PolicyEffectOnZ[1]
        b3_3pct<-PolicyEffectOnZ[2]
        b4zz_1s<-PolicyEffectOnZ[3]
        b4zz_1o<-PolicyEffectOnZ[4]
        #b4x_5pct<-PolicyEffectOnZ[3]
        #policy on X
        #PolicyEffectOnX<-c(0.015,0.03,0.05)
        #PolicyEffectOnX<-c(log(1.05),log(1.10),log(1.15)) 
        b1_1Pfivepct<-PolicyEffectOnX[1]
        b2_1P3pct<-PolicyEffectOnX[2]
        b3_1P5pct<-PolicyEffectOnX[3]
        #X_parameter
        #paraneterX<-c(log(2000),log(1),log(0.98)) #Bo, trendB1,B effect on X
        BO_x2<-paraneterX[1]
        B1_x2<-paraneterX[2]
        B_x2<-paraneterX[3]
        
        
        
        ##varyingZpretrend
        #DifferentZSlope<-c(-0.0035,-0.0020,-0.0005)
        #if (zslopPercent<=0){ ZmStrongerparaV_zfaster StrongerparalleZ_Vp_zfaster
        mildparalleZ_Vf<-DifferentZSlope[1]
        StrongparalleZ_Vf<-DifferentZSlope[2]
        StrongerparalleZ_Vf<-DifferentZSlope[3]
        StrongerparalleZ_VfInc<-DifferentZSlope[4]
        #}else {
        mildparalleZ_Vp<-((100-zslopPercent)/100)*B1
        StrongparalleZ_Vp<-((100-(zslopPercent+PincZsl))/100)*B1
        StrongerparalleZ_Vp<-((100-(zslopPercent+(2*PincZsl)))/100)*B1
        StrongerparalleZ_Vp_zfaster<-((100+(zslopPercent+(2*PincZsl)))/100)*B1
        #}#Print the slope being used
        
        if (zslopPercent<=0){
          #print(mildparalleZ_Vf)
          #print(StrongparalleZ_Vf)
          #print(StrongerparalleZ_Vf)
          #print(StrongerparalleZ_VfInc)
        }else{
          #print(mildparalleZ_Vp)
          #print(StrongparalleZ_Vp)
          #print(StrongerparalleZ_Vp)
          #print(StrongerparalleZ_Vp_zfaster)
        }
        
        
        
        T<-n #100 for now
        t  <- 1:T
        tquad<-t^2
        #t0 <- 25
        if (T%%2==0){
          t0<-((T/2)+1)
          
        }else{
          t0<- floor((T/2)+0.5) 
        }
        
        #print(t0)
        
        P  <- as.integer(t >= t0)     # policy step at t=25
        #N  <- rep(1000, T)            # exposure (can vary if you want)
        
        ##letting confounder have level drop 3 step early around policy
        tq<-t0-lagtq
        #print(tq)
        
        #seasonalX<-c(log(1.04),log(0.98),log(1.04),log(1.01)) #Just picked berneal et al to inform this seasonality component
        Season<-seasonalX[1]*sin((2*pi*t*1)/12)+seasonalX[2]*cos((2*pi*t*1)/12)+seasonalX[3]*sin((2*pi*t*2)/12)+seasonalX[4]*cos((2*pi*t*2)/12)
        
        
        
        ## ----- CONFOUNDERS (plausible, not collinear with policy) -----
        # X1: near-coincident shock at t>=22 + small drift + noise lagpstToX1
        #X1 <- x1drop_schock* as.integer(t >=tq & (t<=(t0+lagtq))) + driftX1* (t)+Season + rnorm(T, 0, 0.20)
        X1 <- x1drop_schock* as.integer(t >=tq & (t<=(t0+lagpstT0X1))) + driftX1* (t)+Season + rnorm(T, 0, 0.20)
        X1 <- X1 - mean(X1[t < t0])
        #X1<- as.integer(t >= (t0 - 3))
        # Seasonality (optional but realistic)
        #S1 <- sin(2*pi*t/12)
        #C1 <- cos(2*pi*t/12)
        ## --- log-mean (eta) and mean (mu) ---
        #X2 count
        #policy on X
        
        #PolicyEffectOnX<-c(log(1.05),log(1.10),log(1.15)) #let me try 5% 10% 15% increase on X
        ##X parameters
        #paraneterX<-c(log(2000),log(1.001),log(0.98))
        #PolicyEffectOnX<-c(0.015,0.03,0.05)
        #b1_1Pfivepct<-PolicyEffectOnX[1]
        #b2_1P3pct<-PolicyEffectOnX[2]
        #b3_1P5pct<-PolicyEffectOnX[3]
        #etaX2_1pointFivepct<-BO_x2+ B1_x2*t +b1_1Pfivepct*P #remove time
        #paraneterX<-c(log(2000),log(1),log(0.98)) #Bo, trendB1,B effect on X
        #BO_x2<-paraneterX[1]
        #B1_x2<-paraneterX[2]
        #B_x2<-paraneterX[3]
        #X2 soem of effect is throught indirect effect
        #.__rng_before <- .Random.seed
        #seed_before <- .Random.seed
        if (genX2==TRUE) with_preserve_seed({
          #seed_before <- .Random.seed
          #set.seed(seed + i)
          T_total   <- B2                      # Interpret vcoeffsY[jek] as TOTAL
          p_share   <- prop_MedX2
          #beta_X2   <- log(0.98)               # −2% per +1 scaled unit (per +1000 visits)
          IE        <- p_share * T_total       # target indirect on log scale
          s         <- (IE /B_x2) * (tensX2 / meanpreX2)   # raw % step on X2 due to policy
          DE        <- T_total - IE            # direct policy effect on Y
          
          # --- Generate raw X2 with a post-policy step of size 's' ---
          etaX2 <- BO_x2 + B1_x2 * t + log(1 + s) * as.integer(t >= (t0 + lagX2)) + rnorm(T, 0, 0.20)
          muX2  <- exp(etaX2)
          X2_raw <- rpois(T, lambda = muX2)
          
          # Pre-policy mean (you can keep using meanpreX2 if you prefer it fixed)
          mu_pre_obs <- mean(X2_raw[t < t0])
          
          # Scaled mediator used in Y model (center + divide by tensX2)
          X2_scaled <- (X2_raw - mu_pre_obs) / tensX2
          #.Random.seed <- seed_before
        })else{
          #print(" ")
        }
        #.Random.seed <- seed_before
        
        etaX2_1pointFivepct<-BO_x2+B1_x2*t+b1_1Pfivepct*(as.integer(t>=(t0+lagX2)))+rnorm(T, 0, 0.20)
        mu  <- exp(etaX2_1pointFivepct)
        X2_1PointFivepct<- rpois(T, mu)
        x2RawX2_1PointFivepct<-X2_1PointFivepct
        X2_1PointFivepct<- (X2_1PointFivepct- mean(X2_1PointFivepct[t < t0])) / sd(X2_1PointFivepct[t < t0])
        
        
        #etaX2_3pct<-BO_x2+ B1_x2*t +b2_1P3pct*P
        etaX2_3pct<-BO_x2+B1_x2*t+b2_1P3pct *(as.integer(t>=(t0+lagX2)))+rnorm(T, 0, 0.20)
        mu  <- exp(etaX2_3pct)
        ## --- simulate counts ---
        X2_3pct<- rpois(T, mu)                 # Poisson
        x2RawX2_3pct<-X2_3pct
        X2_3pct<- (X2_3pct- mean(X2_3pct[t < t0])) / sd(X2_3pct[t < t0])
        
        ##etaX2_5pct<-BO_x2+ B1_x2*t +b3_1P5pct*P  
        etaX2_5pct<-BO_x2+B1_x2*t+b3_1P5pct*(as.integer(t>=(t0+lagX2)))+rnorm(T, 0, 0.20)
        mu  <- exp(etaX2_5pct)
        ## --- simulate counts ---
        X2_5pct<- rpois(T, mu)                 # Poisson
        x2rawX2_5pct<-X2_5pct
        X2_5pct<- (X2_5pct- mean(X2_5pct[t < t0])) / sd(X2_5pct[t < t0])
        
        
        
        
        ## Quick checks for (non-)collinearity with policy
        cor_P_X1 <- cor(P, X1)   # should be far from 1
        #print(round(cor_P_X1, 3))
        
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
        #rho    <-rhofd        #0.4 # choose from {0, 0.2, 0.4, 0.6, 0.8}/I took the mean Turner
        sigma2 <- sigma2fd
        #gen_pair_ar1 <- function(T, rho_y, rho_z, sigma2, alpha)
        # --- shared + idiosyncratic AR(1) components ---
       rho_y<-rof
       rho_z_var<-rof
       alpha<-alpha
       rho_z_const<-rho_z_const
       ZY_AR<-gen_y_twoz_ar1 (T,
                                  rho_y,
                                  rho_z_var,
                                  rho_z_const = 0.4,
                                  sigma2,
                                  alpha)
        
        # --- mix to achieve desired shared confounding level ---
        #u_yt <- shared_conf * u_shared + sqrt(1 - shared_conf^2) * u_y_ind
        #u_ct <- shared_conf * u_shared + sqrt(1 - shared_conf^2) * u_z_ind
        
        u_yt<-ZY_AR$u_yt
        u_ct<-ZY_AR$u_ct
        u_ctconstant<-ZY_AR$u_ctconstant
        ##extact innovations
        e_y <- ZY_AR$e_y
        e_z <- ZY_AR$e_z
        #})
        ##AR for control
        ##constant magnitude control series
        
        if (u_s==1){
          ## ----- GENERATE COUNTS (NB2) -----
          # Linear predictors (log-means)
          ##common trend data
          ##Y affected by Z
          #vcoeffsY<-c(log(700),log(0.995),log(0.96),log(0.70),log(0.60),(log(0.84)/(x1drop_schock)),log(1),350)
          
          etaY <- B0 + B1 * t + B2 * P + B3 * X1+B4*(t^2)+u_yt
          muY  <- exp(etaY)
          Y    <- rnbinom(T, size = thetaY, mu = muY)
          
          etaZ <- b0 + b1 * t+ b2* X1+b4*(t^2)+u_ct
          muZ  <- exp(etaZ)
          Z    <- rnbinom(T, size = thetaZ, mu = muZ)
          
          mZ_zeropoint4Cor<-b0 + b1 * t+ b2* X1+b4*(t^2)+u_ctconstant
          muz<-exp(mZ_zeropoint4Cor)
          Z_zeropoint4Cor<-rnbinom(T, size = thetaZ, mu =muz)
          ##Cofficient Z
          #b0<-vcoeffsZ[1]
          #b1<-vcoeffsZ[2]
          #b2<-vcoeffsZ[3]
          #b3<-vcoeffsZ[4]
          #b4<-vcoeffsZ[5] #carvature
          #thetaZ<-vcoeffsZ[6]
          #vcoeffsZ<-c(log(900),log(0.995),(log(0.84)/(x1drop_schock)),log(1.05),log(1),350) # bo,b1(time),b2(confounder),thetaZ.Lets tKE 90% trend,coef
          #Spillover:
          
          ##Y with X2 that is affected by policy
          # after simulating X2
          #xbar_pre <- mean(X2[t < t0])
          #X2c      <- X2 - xbar_pre          # raw units, just shifted
          #beta_X2  <-B_x2/ 100        # 2% drop per +100 visits
          beta_X2  <-B_x2 #X2_centred
          
          etaY <- B0 + B1*t + B2*P + B3*X1 + B4*(t^2) + beta_X2*X2_1PointFivepct + u_yt  #Policy increases X by 2%
          muY  <- exp(etaY)
          YPolicyon_X1pFivepct<- rnbinom(T, size = thetaY, mu = muY)
          
          
          
          etaY <- B0 + B1*t + B2*P + B3*X1 + B4*(t^2) + beta_X2*X2_3pct+ u_yt  #Policy increases X by 2%
          muY  <- exp(etaY)
          YPolicyon_X3pct<- rnbinom(T, size = thetaY, mu = muY)
          
          
          etaY <- B0 + B1*t + B2*P + B3*X1 + B4*(t^2) + beta_X2*X2_5pct+ u_yt  #Policy increases X by 2%
          muY  <- exp(etaY)
          YPolicyon_X5pct<- rnbinom(T, size = thetaY, mu = muY) 
          
          
          
          
          #SPILL OVER
          ##Spillover effect(opposite direction)
          etaZ <- b0 + b1 * t+b2* X1+((b4zz_1o)*P)+b4*(t^2)+u_ct
          muZ  <- exp(etaZ)
          ZopposieDirSpil<- rnbinom(T, size = thetaZ, mu = muZ)
          #same direction
          etaZ <- b0 + b1 * t+b2* X1+((b4zz_1s)*P)+b4*(t^2)+u_ct   
          muZ  <- exp(etaZ)
          ZsameDirSpil<- rnbinom(T, size = thetaZ, mu = muZ)
          #b3_1Pfivepct<-PolicyEffectOnZ[1]
          #b3_3pct<-PolicyEffectOnZ[2] b3_3pct
          
          
          etaZ <- b0 + b1 * t+b2* X1+(b3_1Pfivepct*(P))+b4*(t^2)+u_ct   
          muZ  <- exp(etaZ)
          ZsameDirSpil_1Pfivepct<- rnbinom(T, size = thetaZ, mu = muZ)
          
          etaZ <- b0 + b1 * t+b2* X1+(b3_3pct*(P))+b4*(t^2)+u_ct   
          muZ  <- exp(etaZ)
          ZsameDirSpil_3pct<- rnbinom(T, size = thetaZ, mu = muZ)
          
          ##TREND VIOLATION PLAYING WITH TIME(Z declining slower)
          if (zslopPercent<=0){
            etaZ <- b0 + mildparalleZ_Vf * t+ b2* X1+b4*(t^2)+u_ct
            muZ  <- exp(etaZ)
            ZmildparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongparalleZ_Vf * t+ b2* X1+b4*(t^2)+u_ct
            muZ  <- exp(etaZ)
            ZmStrongparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongerparalleZ_Vf* t+ b2* X1+b4*(t^2)+u_ct
            muZ  <- exp(etaZ)
            ZmStrongerparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            
            with_preserve_seed({
              etaZ <- b0 +StrongerparalleZ_VfInc* t+ b2* X1+b4*(t^2)+u_ct
              muZ  <- exp(etaZ)
              ZmStrongerparaV_zfaster<- rnbinom(T, size = thetaZ, mu = muZ)
            }
            
            )
            
          }
          else{
            with_preserve_seed({
              etaZ <- b0 + mildparalleZ_Vp * t+ b2* X1+b4*(t^2)+u_ct
              muZ  <- exp(etaZ)
              ZmildparaV<- rnbinom(T, size = thetaZ, mu = muZ)
              
              etaZ <- b0 + StrongparalleZ_Vp * t+ b2* X1+b4*(t^2)+u_ct
              muZ  <- exp(etaZ)
              ZmStrongparaV<- rnbinom(T, size = thetaZ, mu = muZ)
              
              etaZ <- b0 + StrongerparalleZ_Vp* t+ b2* X1+b4*(t^2)+u_ct
              muZ  <- exp(etaZ)
              ZmStrongerparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            })
            
            with_preserve_seed({
              etaZ <- b0 +StrongerparalleZ_Vp_zfaster* t+ b2* X1+b4*(t^2)+u_ct
              muZ  <- exp(etaZ)
              ZmStrongerparaV_zfaster<- rnbinom(T, size = thetaZ, mu = muZ)
            }
            
            )
            
            etaZ <- b0 + mildparalleZ_Vf * t+ b2* X1+b4*(t^2)+u_ct
            muZ  <- exp(etaZ)
            rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongparalleZ_Vf * t+ b2* X1+b4*(t^2)+u_ct
            muZ  <- exp(etaZ)
            rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongerparalleZ_Vf* t+ b2* X1+b4*(t^2)+u_ct
            muZ  <- exp(etaZ)
            rnbinom(T, size = thetaZ, mu = muZ)
            
            
          }
          
          ##unshared Z confounder
          #etaZ <- b0 + b1 * t+b4*(t^2)+u
          #muZ  <- exp(etaZ)
          #Z_unXrm<- rnbinom(T, size = thetaZ, mu = muZ)
          
          #etaZ <- b0 + b1 * t+b2* X1+u
          #muZ  <- exp(etaZ)
          #Z_unCurve<- rnbinom(T, size = thetaZ, mu = muZ)
          ##Y mediated
          if (genX2==TRUE) with_preserve_seed({
            #seed_before <- .Random.seed
            #set.seed(seed + i)
            etaYM<-B0 + B1*t + DE*P + B3*X1 + B4*(t^2) + beta_X2 * X2_scaled+u_yt #Policy increases X by 2%
            muY<- exp(etaYM)
            Y_aff<- rnbinom(T, size = thetaY, mu = muY)
            dat<- data.frame(t,Z_zeropoint4Cor,# >>> ADD THESE <<<
                             u_yt = u_yt,
                             u_ct = u_ct,
                             e_y  = e_y,
                             e_z  = e_z,ZmStrongerparaV_zfaster,Y_aff,X2_scaled,x2RawX2_3pct,x2rawX2_5pct,x2RawX2_1PointFivepct, P,X1,Y,YPolicyon_X1pFivepct,YPolicyon_X3pct,YPolicyon_X5pct,X2_1PointFivepct,X2_3pct,X2_5pct,ZmildparaV,ZmStrongparaV,ZmStrongerparaV,ZsameDirSpil_1Pfivepct,ZsameDirSpil_3pct,Z,ZsameDirSpil,ZopposieDirSpil,t0=t0,tq=tq)
            #.Random.seed <- seed_before
          })else{
            dat<- data.frame(t,Z_zeropoint4Cor,# >>> ADD THESE <<<
                             u_yt = u_yt,
                             u_ct = u_ct,
                             e_y  = e_y,
                             e_z  = e_z,ZmStrongerparaV_zfaster,x2RawX2_3pct,x2rawX2_5pct,x2RawX2_1PointFivepct, P,X1,Y,YPolicyon_X1pFivepct,YPolicyon_X3pct,YPolicyon_X5pct,X2_1PointFivepct,X2_3pct,X2_5pct,ZmildparaV,ZmStrongparaV,ZmStrongerparaV,ZsameDirSpil_1Pfivepct,ZsameDirSpil_3pct,Z,ZsameDirSpil,ZopposieDirSpil,t0=t0,tq=tq)
          }
          
          
          
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
          #etaY <- B0 + B1 * t + B2 * P + B3 * X1+B4*(t^2)+u
          etaY <- B0 + B1 * t + B2 * P + B3 * X1+B4*(t^2)
          muY  <- exp(etaY)
          Y    <- rnbinom(T, size = thetaY, mu = muY)
          
          #etaZ <- b0 + b1 * t+ b2* X1+b4*(t^2)+u
          etaZ <- b0 + b1 * t+ b2* X1+b4*(t^2)
          muZ  <- exp(etaZ)
          Z    <- rnbinom(T, size = thetaZ, mu = muZ)
          
          
          mZ_zeropoint4Cor<-b0 + b1 * t+ b2* X1+b4*(t^2)
          muz<-exp(mZ_zeropoint4Cor)
          Z_zeropoint4Cor<-rnbinom(T, size = thetaZ, mu =muz)
          ##Y with X2 that is affected by policy
          # after simulating X2
          #xbar_pre <- mean(X2[t < t0])
          #X2c      <- X2 - xbar_pre          # raw units, just shifted
          beta_X2  <-B_x2       # 2% drop per +100 visits
          
          etaY <- B0 + B1*t + B2*P + B3*X1 + B4*(t^2) + beta_X2*X2_1PointFivepct  #Policy increases X by 2%
          muY  <- exp(etaY)
          YPolicyon_X1pFivepct<- rnbinom(T, size = thetaY, mu = muY)
          
          
          
          etaY <- B0 + B1*t + B2*P + B3*X1 + B4*(t^2) + beta_X2*X2_3pct #Policy increases X by 2%
          muY  <- exp(etaY)
          YPolicyon_X3pct<- rnbinom(T, size = thetaY, mu = muY)
          
          
          etaY <- B0 + B1*t + B2*P + B3*X1 + B4*(t^2) + beta_X2*X2_5pct  #Policy increases X by 2%
          muY  <- exp(etaY)
          YPolicyon_X5pct<- rnbinom(T, size = thetaY, mu = muY)
          
          #SPILL OVER
          ##Spillover effect(opposite direction)
          #etaZ <- b0 + b1 * t+b2* X1+b3*P+b4*(t^2)
          #opposite
          #PolicyEffectOnZ<-c(-0.015,-0.03,-0.05)
          #PolicyEffectOnZ<-c(log(0.95),log(0.90),log(0.85),log(1.15))
          #b3_1Pfivepct<-PolicyEffectOnZ[1]
          #b3_3pct<-PolicyEffectOnZ[2]
          #b4zz_1s<-PolicyEffectOnZ[3]
          #b4zz_1o<-PolicyEffectOnZ[4]
          etaZ <- b0 + b1 * t+b2* X1+((b4zz_1o)*P)+b4*(t^2)
          muZ  <- exp(etaZ)
          ZopposieDirSpil<- rnbinom(T, size = thetaZ, mu = muZ)
          #same direction
          etaZ <- b0 + b1 * t+b2* X1+((b4zz_1s)*P)+b4*(t^2) 
          muZ  <- exp(etaZ)
          ZsameDirSpil<- rnbinom(T, size = thetaZ, mu = muZ)
          
          
          
          etaZ <- b0 + b1 * t+b2* X1+(b3_1Pfivepct*(P))+b4*(t^2) 
          muZ  <- exp(etaZ)
          ZsameDirSpil_1Pfivepct<- rnbinom(T, size = thetaZ, mu = muZ)
          
          etaZ <- b0 + b1 * t+b2* X1+(b3_3pct*(P))+b4*(t^2)  
          muZ  <- exp(etaZ)
          ZsameDirSpil_3pct<- rnbinom(T, size = thetaZ, mu = muZ)
          
          
          ##TREND VIOLATION PLAYING WITH TIME(Z declining slower)
          if (zslopPercent<=0){
            etaZ <- b0 + mildparalleZ_Vf * t+ b2* X1+b4*(t^2)
            muZ  <- exp(etaZ)
            ZmildparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongparalleZ_Vf * t+ b2* X1+b4*(t^2)
            muZ  <- exp(etaZ)
            ZmStrongparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongerparalleZ_Vf* t+ b2* X1+b4*(t^2)
            muZ  <- exp(etaZ)
            ZmStrongerparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            
            with_preserve_seed({
              etaZ <- b0 +StrongerparalleZ_VfInc* t+ b2* X1+b4*(t^2)
              muZ  <- exp(etaZ)
              ZmStrongerparaV_zfaster<- rnbinom(T, size = thetaZ, mu = muZ)
            }
            
            )
            
          }
          else{
            with_preserve_seed({
              etaZ <- b0 + mildparalleZ_Vp * t+ b2* X1+b4*(t^2)
              muZ  <- exp(etaZ)
              ZmildparaV<- rnbinom(T, size = thetaZ, mu = muZ)
              
              etaZ <- b0 + StrongparalleZ_Vp * t+ b2* X1+b4*(t^2)
              muZ  <- exp(etaZ)
              ZmStrongparaV<- rnbinom(T, size = thetaZ, mu = muZ)
              
              etaZ <- b0 + StrongerparalleZ_Vp* t+ b2* X1+b4*(t^2)
              muZ  <- exp(etaZ)
              ZmStrongerparaV<- rnbinom(T, size = thetaZ, mu = muZ)
            })
            
            with_preserve_seed({
              etaZ <- b0 +StrongerparalleZ_Vp_zfaster* t+ b2* X1+b4*(t^2)
              muZ  <- exp(etaZ)
              ZmStrongerparaV_zfaster<- rnbinom(T, size = thetaZ, mu = muZ)
            }
            
            )
            
            etaZ <- b0 + mildparalleZ_Vf * t+ b2* X1+b4*(t^2)
            muZ  <- exp(etaZ)
            rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongparalleZ_Vf * t+ b2* X1+b4*(t^2)
            muZ  <- exp(etaZ)
            rnbinom(T, size = thetaZ, mu = muZ)
            
            etaZ <- b0 + StrongerparalleZ_Vf* t+ b2* X1+b4*(t^2)
            muZ  <- exp(etaZ)
            rnbinom(T, size = thetaZ, mu = muZ)
            
            
          }
          
          
          ##unshared Z confounder
          #etaZ <- b0 + b1 * t+b4*(t^2)
          #muZ  <- exp(etaZ)
          #Z_unXrm<- rnbinom(T, size = thetaZ, mu = muZ)
          
          #etaZ <- b0 + b1 * t+b2* X1
          #muZ  <- exp(etaZ)
          #Z_unCurve<- rnbinom(T, size = thetaZ, mu = muZ)
          ##Y mediated
          if (genX2==TRUE) with_preserve_seed({
            #seed_before <- .Random.seed
            #set.seed(seed + i)
            etaYM<-B0 + B1*t + DE*P + B3*X1 + B4*(t^2) + beta_X2 * X2_scaled #Policy increases X by 2%
            muY<- exp(etaYM)
            Y_aff<- rnbinom(T, size = thetaY, mu = muY)
            dat<- data.frame(t,Z_zeropoint4Cor,ZmStrongerparaV_zfaster,Y_aff,X2_scaled,x2RawX2_3pct,x2rawX2_5pct,x2RawX2_1PointFivepct, P,X1,Y,YPolicyon_X1pFivepct,YPolicyon_X3pct,YPolicyon_X5pct,X2_1PointFivepct,X2_3pct,X2_5pct,ZmildparaV,ZmStrongparaV,ZmStrongerparaV,ZsameDirSpil_1Pfivepct,ZsameDirSpil_3pct,Z,ZsameDirSpil,ZopposieDirSpil,t0=t0,tq=tq)
            #.Random.seed <- seed_before
          })else{
            dat<- data.frame(t,Z_zeropoint4Cor,ZmStrongerparaV_zfaster,x2RawX2_3pct,x2rawX2_5pct,x2RawX2_1PointFivepct, P,X1,Y,YPolicyon_X1pFivepct,YPolicyon_X3pct,YPolicyon_X5pct,X2_1PointFivepct,X2_3pct,X2_5pct,ZmildparaV,ZmStrongparaV,ZmStrongerparaV,ZsameDirSpil_1Pfivepct,ZsameDirSpil_3pct,Z,ZsameDirSpil,ZopposieDirSpil,t0=t0,tq=tq)
          }
          
          
          
          #dat<- data.frame(t,Y_aff,X2_scaled ,x2RawX2_3pct,x2rawX2_5pct,x2RawX2_1PointFivepct, P,X1,Y,YPolicyon_X1pFivepct,YPolicyon_X3pct,YPolicyon_X5pct,X2_1PointFivepct,X2_3pct,X2_5pct,ZmildparaV,ZmStrongparaV,ZmStrongerparaV,ZsameDirSpil_1Pfivepct,ZsameDirSpil_3pct,Z,ZsameDirSpil,ZopposieDirSpil,t0=t0,tq=tq)
          dat$j<-i
          sim_listcommon[[i]]<-dat
          ##unshared trend data
          #etaY <- B0 + B1 * t + B2 * P + g1 * X1+u
          #muY  <- exp(etaY)
          #Y    <- rnbinom(T, size = thetaY, mu = muY)
        }
        
      }
      ## return BOTH
      #list(common = sim_listcommon, uncommon = sim_listUncommon)
      #sim_listcommonBO_Expected<-list()
      #sim_listcommonBO_MinEf<-list()
      #sim_listcommonBO_MaxEf<-list() 
      if (oyaoya==1){
        sim_listcommonBO_MinEf=sim_listcommon
        #print("Minimal effect have",nrow())
      }else if(oyaoya==2){
        sim_listcommonBO_Expected=sim_listcommon 
      }else{
        sim_listcommonBO_MaxEf=sim_listcommon 
      }
      
      
      
      oyaoya=oyaoya+1  
    }
    # ---- ADD: tag each data frame with rho and n (minimal, base R only) ----
    # ---- ADD: tag each data frame with rho and n (robust to transform() scoping) ----
    #cat("corelation is", rof, "\n")
    #cat("Length of series is", jojog, "\n")
    
    rho_tag <- as.character(rof)
    n_tag   <- as.integer(jojog)
    
    add_tags <- function(d) {
      d$rho <- rho_tag
      d$n   <- n_tag
      d
    }
    
    sim_listcommonBO_MinEf    <- lapply(sim_listcommonBO_MinEf,    add_tags)
    sim_listcommonBO_Expected <- lapply(sim_listcommonBO_Expected, add_tags)
    sim_listcommonBO_MaxEf    <- lapply(sim_listcommonBO_MaxEf,    add_tags)
    
    res[[paste0(as.character(rof), "_", jojog)]] <- list(
      sim_listcommonBO_MinEf    = sim_listcommonBO_MinEf,
      sim_listcommonBO_Expected = sim_listcommonBO_Expected,
      sim_listcommonBO_MaxEf    = sim_listcommonBO_MaxEf
    )
    
    #list(sim_listcommonBO_MinEf=sim_listcommonBO_MinEf, sim_listcommonBO_Expected=sim_listcommonBO_Expected,sim_listcommonBO_MaxEf=sim_listcommonBO_MaxEf)
    ##Print rho and n here rho jojog
    #garbe collection
    rm(sim_listcommonBO_MinEf,sim_listcommonBO_Expected,sim_listcommonBO_MaxEf)
    tmp_gc <- gc()
  print(tmp_gc)
  rm(tmp_gc)
  }
  }
return(res)

  ##end function
  }
  
##running this function to generate shared and unshared confounder# run u_s=0 for not adding AR(1)
# Run with identical seed inside the function for fairness
##COMAPRE SIMULATED Y in 10 simulations(nsim,seasonalX,vcoeffsY,vcoeffsZ,rhofd, sigma2fd,n,u_s=0,genX2=TRUE,zslopPercent=5,PincZsl=5,prop_MedX2=0.2,tensX2=1000,meanpreX2=2000,seed=123){ ##Let var be 3 times the mean. mean to variance ration for choosing overdispersion.Chose mean to var ratio of 2;mean=B0
#setting seed
library(dplyr)

## --- Set the pair you want to test ---
rho_key <- "0.4"    # choose "0", "0.2", "0.4", "0.6", "0.8"
n_key   <- "150"    # choose the n you looped (e.g., "100", "150", ...)
#rho_key <- "0.4"
#n_key   <- "150"
## 1) RUN with genX2 = FALSE
o_false <- fredSiProject(
  10, seasonalX, vcoeffsY, vcoeffsZ,
  rhofd, Rhosingle = 0.4, SingleRho = TRUE,
  0.1, as.integer(n_key), u_s = 1, genX2 = FALSE
)

## Pull the correct object created by your assign():
#obj_name <- paste0("sim_listcommonBO_", rho_key, "_", n_key)
#f0_list  <- get(obj_name)$sim_listcommonBO_Expected

obj_key <- paste0(rho_key, "_", n_key)   # e.g. "0.4_150"

f0_list <- o_false[[obj_key]]$sim_listcommonBO_Expected


## Bind and keep only what we’ll compare
f0 <- dplyr::bind_rows(f0_list, .id = "sim") |>
  dplyr::arrange(as.integer(sim), t) |>
  dplyr::select(sim, j, t, Y)

## 2) RUN with genX2 = TRUE  (this will overwrite the same global name)
o_true <- fredSiProject(
  10, seasonalX, vcoeffsY, vcoeffsZ,
  rhofd, Rhosingle = 0.4, SingleRho = TRUE,
  0.1, as.integer(n_key), u_s = 1, genX2 = TRUE
)

## Grab the same rho/n block again (now from the genX2=TRUE run)
obj_key <- paste0(rho_key, "_", n_key)   # e.g. "0.4_150"

f1_list <-o_true[[obj_key]]$sim_listcommonBO_Expected

#f1_list <- get(obj_name)$sim_listcommonBO_Expected
f1 <- dplyr::bind_rows(f1_list, .id = "sim") |>
  dplyr::arrange(as.integer(sim), t) |>
  dplyr::select(sim, j, t, Y)

## 3) Compare Y
identical(f0$Y, f1$Y)


rm(f0,f1,f0_list,o_false,o_true ,f1_list )
gc()

##Check if Y and Z dont change
## --- choose the block you want to compare ---
rho_key <- "0.4"
n_key   <- "150"
obj_name <- paste0(rho_key, "_", n_key)  # <-- use the combined key

## Run A (percent violation on Z slopes)
outA <- fredSiProject(
  10, seasonalX, vcoeffsY, vcoeffsZ,
  rhofd, Rhosingle = 0.4, SingleRho = TRUE,
  0.1, 150, u_s = 1, genX2 = TRUE,
  zslopPercent = 5, PincZsl = 5, averaging_n = 150
)

## Run B (fixed Z slopes)
outB <- fredSiProject(
  10, seasonalX, vcoeffsY, vcoeffsZ,
  rhofd, Rhosingle = 0.4, SingleRho = TRUE,
  0.1, 150, u_s = 1, genX2 = TRUE,
  zslopPercent = 0, averaging_n = 150
)

## ---- Baselines should be identical (Y, Z) ----
A <- do.call(rbind, outA[[obj_name]]$sim_listcommonBO_MinEf)[, c("j","t","Y","Z")]
B <- do.call(rbind, outB[[obj_name]]$sim_listcommonBO_MinEf)[, c("j","t","Y","Z")]

A <- A[order(A$j, A$t), ]
B <- B[order(B$j, B$t), ]

identical(A$Y, B$Y)  # expect TRUE
identical(A$Z, B$Z)  # expect TRUE

## ---- Trend-violation series should differ ----
VA <- do.call(rbind, outA[[obj_name]]$sim_listcommonBO_Expected)[, c("j","t","ZmildparaV")]
VB <- do.call(rbind, outB[[obj_name]]$sim_listcommonBO_Expected)[, c("j","t","ZmildparaV")]

VA <- VA[order(VA$j, VA$t), ]
VB <- VB[order(VB$j, VB$t), ]

identical(VA$ZmildparaV, VB$ZmildparaV)  # expect FALSE



rm(outA,outB,A,B,VA,VB)
gc()
##remove
#rm(list = ls())
##RUN MAIN simulation(nsim,seasonalX,vcoeffsY,vcoeffsZ,rhofd, sigma2fd,n,u_s=0,genX2=TRUE,zslopPercent=5,PincZsl=5,prop_MedX2=0.2,tensX2=1000,meanpreX2=2000,seed=123){ ##Let var be 3 times the mean. mean to variance ration for choosing overdispersion.Chose mean to var ratio of 2;mean=B0
#setting seed
#averaging_n=150,you can set averaging n to any values for which you want to preent results if results are so many:I will do 7000 but let me test with 200
#8688
out <-fredSiProject(20,seasonalX,vcoeffsY,vcoeffsZ,rhofd,Rhosingle=0.4,SingleRho=FALSE,0.1,n,u_s =1,genX2=TRUE,zslopPercent=10,PincZsl=10,averaging_n=150,shared_conf =0.85)
names(out)
#out <-fredSiProject(7000,seasonalX,vcoeffsY,vcoeffsZ,0.4,0.1,150,u_s =0)
#names(out)
#d1_common   <- out$common[[1]]
#d1_uncommon <- out$uncommon[[1]]
##Bind all shared confounder simulation
#library(dplyr)

#rho_key <- "0.4"  # choose the rho you want
n_key <- "150" 

for (eve in rhofd ){
  #foget<-as.character(eve)
  if (eve==0.4){
    n_key <- "150" 
    rho_key<-as.character(eve)
    obj_name<- paste0(rho_key, "_", n_key)
    # 1) Grab ONLY the Expected sims for rho = 0.4
    expected_list <- out[[obj_name]]$sim_listcommonBO_Expected
    
    # 2) Bind all sims, keep an id if you want (drop `.id="sim"` if not needed)
    common_all <- dplyr::bind_rows(expected_list, .id = "sim") %>%
      arrange(as.integer(sim), j, t)
    
    # 3) Pull t0/tq from the first row (they’re constant within a run)
    t0 <- common_all$t0[1]
    tq <- common_all$tq[1]
    print(t0); print(tq)
    
    # 4) Add derived columns and Estimand flag
    common_all <- common_all %>%
      mutate(
        tce = t - t0,
        Policy_Time = tce * P,
        EstimandScenario =1
      )
    
    # 5) Save
    common_all$sim<-NULL
    if (!dir.exists("data")) dir.create("data", recursive = TRUE)
    #saveRDS(common_all, "data/common_allExp.rds")
    
    
    
    
    # 1) Grab ONLY the Expected sims for rho = 0.4
    expected_list <- out[[obj_name]]$sim_listcommonBO_MinEf
    
    # 2) Bind all sims, keep an id if you want (drop `.id="sim"` if not needed)
    common_allMin<- dplyr::bind_rows(expected_list, .id = "sim") %>%
      arrange(as.integer(sim), j, t)
    
    # 3) Pull t0/tq from the first row (they’re constant within a run)
    t0 <- common_allMin$t0[1]
    tq <- common_allMin$tq[1]
    print(t0); print(tq)
    
    # 4) Add derived columns and Estimand flag
    common_allMin<- common_allMin %>%
      mutate(
        tce = t - t0,
        Policy_Time = tce * P,
        EstimandScenario =2
      )
    common_allMin$sim<-NULL
    # 5) Save
    if (!dir.exists("data")) dir.create("data", recursive = TRUE)
    #saveRDS(common_allMin, "data/common_allMin.rds")
    
    
    
    
    # 1) Grab ONLY the Expected sims for rho = 0.4
    expected_list <- out[[obj_name]]$sim_listcommonBO_MaxEf
    
    # 2) Bind all sims, keep an id if you want (drop `.id="sim"` if not needed)
    common_allMax <- dplyr::bind_rows(expected_list, .id = "sim") %>%
      arrange(as.integer(sim), j, t)
    
    # 3) Pull t0/tq from the first row (they’re constant within a run)
    t0 <- common_allMax$t0[1]
    tq <- common_allMax$tq[1]
    print(t0); print(tq)
    
    # 4) Add derived columns and Estimand flag
    common_allMax<- common_allMax %>%
      mutate(
        tce = t - t0,
        Policy_Time = tce * P,
        EstimandScenario =3
      )
    
    # 5) Save
    common_allMax$sim<-NULL
    if (!dir.exists("data")) dir.create("data", recursive = TRUE)
    #saveRDS(common_allMax, "data/common_allMax.rds")
    
    
    
    
    
    #View(common_all)
    ##Append_ThemAll
    SimulatedData<-rbind(common_all,common_allMin, common_allMax)
    SimulatedData$EstimandScenario<-factor(SimulatedData$EstimandScenario,
                                           levels = c(1, 2, 3),
                                           labels = c("Expected", "minimum", "Maximum"))
    #Save dataset.
    saveRDS(SimulatedData,"data/SimulatedData.rds")
    ##END OF SIMULATION
    
    rm(SimulatedData,common_all,common_allMin,common_allMax,expected_list)
    gc()
    
  }
  rho_key<-as.character(eve)
  for (n_key in n) {
    obj_name <- paste0(rho_key, "_", n_key)
    ##Other combinations
    # 1) Grab ONLY the Expected sims for rho = 0.4
    expected_list <- out[[obj_name]]$sim_listcommonBO_Expected
    
    # 2) Bind all sims, keep an id if you want (drop `.id="sim"` if not needed)
    common_all <- dplyr::bind_rows(expected_list, .id = "sim") %>%
      arrange(as.integer(sim), j, t)
    
    # 3) Pull t0/tq from the first row (they’re constant within a run)
    t0 <- common_all$t0[1]
    tq <- common_all$tq[1]
    print(t0); print(tq)
    
    # 4) Add derived columns and Estimand flag
    common_all <- common_all %>%
      mutate(
        tce = t - t0,
        Policy_Time = tce * P,
        EstimandScenario =2
      )
    
    # 5) Save
    common_all$sim<-NULL
    if (!dir.exists("data")) dir.create("data", recursive = TRUE)
    #saveRDS(common_all, "data/common_allExp.rds")
    
    
    
    
    # 1) Grab ONLY the Expected sims for rho = 0.4
    expected_list <- out[[obj_name]]$sim_listcommonBO_MinEf
    
    # 2) Bind all sims, keep an id if you want (drop `.id="sim"` if not needed)
    common_allMin<- dplyr::bind_rows(expected_list, .id = "sim") %>%
      arrange(as.integer(sim), j, t)
    
    # 3) Pull t0/tq from the first row (they’re constant within a run)
    t0 <- common_allMin$t0[1]
    tq <- common_allMin$tq[1]
    print(t0); print(tq)
    
    # 4) Add derived columns and Estimand flag
    common_allMin<- common_allMin %>%
      mutate(
        tce = t - t0,
        Policy_Time = tce * P,
        EstimandScenario =1
      )
    common_allMin$sim<-NULL
    # 5) Save
    if (!dir.exists("data")) dir.create("data", recursive = TRUE)
    #saveRDS(common_allMin, "data/common_allMin.rds")
    
    
    
    
    # 1) Grab ONLY the Expected sims for rho = 0.4
    expected_list <- out[[obj_name]]$sim_listcommonBO_MaxEf
    
    # 2) Bind all sims, keep an id if you want (drop `.id="sim"` if not needed)
    common_allMax <- dplyr::bind_rows(expected_list, .id = "sim") %>%
      arrange(as.integer(sim), j, t)
    
    # 3) Pull t0/tq from the first row (they’re constant within a run)
    t0 <- common_allMax$t0[1]
    tq <- common_allMax$tq[1]
    print(t0); print(tq)
    
    # 4) Add derived columns and Estimand flag
    common_allMax<- common_allMax %>%
      mutate(
        tce = t - t0,
        Policy_Time = tce * P,
        EstimandScenario =3
      )
    
    # 5) Save
    common_allMax$sim<-NULL
    if (!dir.exists("data")) dir.create("data", recursive = TRUE)
    #saveRDS(common_allMax, "data/common_allMax.rds")
    
    
    
    
    
    #View(common_all)
    ##Append_ThemAll
    SimulatedData<-rbind(common_all,common_allMin, common_allMax)
    SimulatedData$EstimandScenario<-factor(SimulatedData$EstimandScenario,
                                           levels = c(1, 2, 3),
                                           labels = c("Small","Moderate","Large"))
    #Save dataset.
    gh<-paste("data/SimulatedData",obj_name,".rds",sep="")
    #SimulatedData$rho<-NULL
    #SimulatedData$rho<-rho_key
    saveRDS(SimulatedData,gh)
    ##END OF SIMULATION
    rm(SimulatedData,common_all,common_allMin,common_allMax,expected_list,gh)
    gc()
    
  }
  
  
} 

##free the memory
rm(out)
gc()


##Checking missing files
rho_keyvec1 <- c("0", "0.2", "0.4", "0.6", "0.8")
n_vec <- c(12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,80,100,150)

missing <- list()

for(rho in rho_keyvec1){
  for(n in n_vec){
    fp <- paste0("data/SimulatedData", rho, "_", n, ".rds")
    if(!file.exists(fp)){
      missing <- c(missing, fp)
    }
  }
}

missing
length(missing)




##Append rho scenarios vs n
## Append rho scenarios vs n
# list of all rho keys
rho_keyvec1 <- c("0", "0.2", "0.4", "0.6", "0.8")

# list of all n values
n_vec <- c(12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,80,100,150)
#     n<-c(12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,80,100,150)
# 1. Read all the RDS files into a list
SimList <- lapply(rho_keyvec1, function(rho) {
  inner_list <- lapply(n_vec, function(n) {
    file_path <- paste0("data/SimulatedData", rho, "_", n, ".rds")
    dat <- readRDS(file_path)
    #dat$rho <- rho     # keep rho
    #dat$n   <- n       # add n
    dat
  })
  dplyr::bind_rows(inner_list)
})

# 2. Combine all into one big data frame
AllRhoCombined <- dplyr::bind_rows(SimList)

# 3. Save the merged dataset
saveRDS(AllRhoCombined, "data/AllRhoCombined_diffAR_ZvsY.rds")
# Export to CSV
#write.csv(AllRhoCombined, "data/AllRhoCombined.csv", row.names = FALSE)
#write_csv(AllRhoCombined, "data/AllRhoCombined.csv")
#write.xlsx(AllRhoCombined, 
           #file = "data/AllRhoCombined.xlsx", 
           #overwrite = TRUE)
# optional: check counts
table(AllRhoCombined$rho,AllRhoCombined$n)
names(AllRhoCombined)

##check corelation between Y and Z I aim for 70%
AllRhoCombined1<-AllRhoCombined[AllRhoCombined$EstimandScenario=="Large" & n==150 & AllRhoCombined$rho==0.8 & AllRhoCombined$j==150 & AllRhoCombined$t<76,]
cor(AllRhoCombined1$Z,AllRhoCombined1$Y)

# Suppose you already created SimulatedData for rho=0.8 and n=150
d_sub <- AllRhoCombined |>
  dplyr::filter(rho == "0.8",
                n == 150,
                EstimandScenario == "Large",
                j == 150,
                t < 76)

# correlations you asked for (latent AR errors + innovations)
cor(d_sub$u_yt, d_sub$u_ct)   # should be close to alpha if rho_y=rho_z=0.8
cor(d_sub$e_y,  d_sub$e_z)    # should be close to alpha

var(d_sub$e_y)                # should be close to sigma2
var(d_sub$e_z)                # should be close to sigma2



rm(AllRhoCombined,SimList)
gc()

##Test some collinearity As you adjust.

common_all<-readRDS("data/SimulatedData.rds")

#common_all<-common_all[common_all$EstimandScenario=="Expected",]
common_all<-common_all[common_all$EstimandScenario=="Expected",]
t0<-common_all$t0[1]
tq<-common_all$tq[1]
print(t0)
print(tq)
##drop these unnecessary columns
common_all$t0<-NULL
common_all$tq<-NULL
nrow(common_all)
#View(common_all)
#Analysis
#delta correction
delta<-0.0001

##ANALYSIS
##plot simulation
## Pre-policy trends overtime.
#Shared confounder common trend loess using first simulated data;
dat<-common_all[common_all$j==1,]
datpre <- subset(dat, P == 0|P == 1)

datpre <- subset(dat, P == 0|P==1)
# 1) All pairwise correlations among the 4 vars
  cor(dat[, c("t", "P", "X1", "X2_1PointFivepct")],
      use = "pairwise.complete.obs", method = "pearson")



##remove
rm(common_all,datpre,dat)
gc()
  

