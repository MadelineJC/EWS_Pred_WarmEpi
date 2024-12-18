require(deSolve) #loading the ODE solver
require(FME)
require(GillespieSSA)


minTemp  <- 273.15  # The minimum temperature for which S_c is calculated
maxTemp  <- 313.15  # The maximum temperature for which S_c is calculated
resol <- 0.01  # Resolution
temp <- seq(from = minTemp, to = maxTemp, by = resol)


# Initialize vectors for each of the temperature-dependent parameters that determine Sc
mu <- vector(mode="numeric", length(temp))
beta <- vector(mode="numeric", length(temp))
alpha <- vector(mode="numeric", length(temp)) 
chi <- vector(mode="numeric", length(temp))
lambda <- vector(mode="numeric", length(temp))
sigma <- vector(mode="numeric", length(temp))
theta <- vector(mode="numeric", length(temp))
omega <- vector(mode="numeric", length(temp))
equil.abundance<- vector(mode="numeric", length(temp))


# Parameters #
# Note when parameterizing that PLOS biol. paper used 15C as T0, and AmNat paper used 20C as T0

# cs is contact allometric scaling, is is is unfection rate allometric scaling

chi0   =  1.255755e+02    ;  sigma0= 5.684252e-09       ;  omega0= 131      ; beta0 <- 2.281060732      ;mu0 <- 0.007507025 # Value at reference temperature                           
E_chi  =  3.073807e-01    ;  E_sigma = 1.106441    ;     E_omega= -0.295  ;  E_beta <- 0.390532437 ;     E_mu <- 0.801482243   # Activation energy (THIS IS WHAT I WANTED TO PLAY AROUND WITH)                          
                            E_sigmaL = 6.766889          ;E_omegaL= 357  ;                              E_muL <- E_mu*5                                                # 'Inactivation' energy at lower threshold (taken as per Molnar et al. 2013)    
E_chiH =  4.196706e+00   ;   E_sigmaH = 5.661312         ;  E_omegaH= 6.32  ; E_betaH <- 1.472373114 ; E_muH <- 4.239515185  # 'Inactivation' energy at upper threshold (taken as per Molnar et al. 2013)          
T0_chi = 293.15 ;      T0_sigma = 293.15        ;  T0_omega= 288.15     ; T0_beta <- 288.15 ; T0_mu <- 288.15  # Reference temperatures (taken as 20C for contact and infection, 15 for mu and omega)                        
                              T_sigmaL = 273.2+14.14539   ;  T_omegaL= 273.2+9.9                  ;   T_muL <-  8.957624193 + 273.15# Lower thresholds                          
T_chiH = 273.15+3.455977e+01  ; T_sigmaH = 273.15+23.43956   ; T_omegaH= 273.15+23.9 ; T_betaH <- 28.8419577 +273.15; T_muH <-  30.23264714 +273.15 #Upper thresholds                        
cs = 0.2574784 ; is = -1.181799;


k <- 8.62*10^(-5)


# Assumptions here #
# 1) Within-host parasite abundance as a function of equilibrium abundance
# Mean infection load across all infected individuals in the experiment was 27.25 clusters #
# Mean equilibrium abundance for the 8 temp treatments is 148.79 clusters #
# Gives proportion of 0.182 on average compared to equilibrium #

abundance.prop = 0.182

# 2) How long the clusters take to burst 
burst.time = 4.5

# 3) How many spores are in each cluster
spores.per.cluster = 12

# 4) How many of the spores are released into the environment from each burst cluster 
burst.size = 0.5

# 5) Daphnia size
size = 2700
mass <-  (0.009*(size*0.001)^2.63)*0.001



# Contact rate is in microliters/min, need to convert to 35L/24 hours to standardize to the volume of the experiment
# There is 3.5x10^7 microliters in 35L, and 1440 minutes in a day.

standardized.contact <- c()

# contact is in microliters / min
for(i in 1:length(temp)){
  standardized.contact[i] <- ((chi0*(mass^cs)*exp((-E_chi/k)*(1/temp[i]-1/T0_chi))/(1+exp(E_chiH/k*(-1/temp[i]+1/T_chiH)))))
  }

# check number on ingestion rate for number of cells in the tank. This is if we fed 350 million in 35L
# Need to calculate ingestion rate (cells/hour), and we have 10,000 algal cells per ml, therefore 10 algal cells per microliter 
ing.cells.rate <- ((standardized.contact)*10)*60

# THIS CONCENTRATION OF ALGAE IS MUCH LOWER THAN IN THE AM NAT PAPER, THEREFORE THE INGESTION OF CELLS IS MUCH LOWER
# AND THE ENSUING GUT RESIDENCE TIME IS WAY LONGER #

# From Kooijman 1993 and Evers and Kooijman: volume of gut in cells = gv (9900) * length^3 (mm^3) / ingested cells rate, this gives gut.res time in hours
gut.res.time <- (9900*((size/1000)^3))/ing.cells.rate


# sigma in this loop is the infection RATE, will calculate sigma probability in step after #


for(t in seq(temp) ){
  chi[t] <- ((((chi0*(mass^cs)*exp((-E_chi/k)*(1/temp[t]-1/T0_chi))/(1+exp(E_chiH/k*(-1/temp[t]+1/T_chiH)))))/3.5E7)*1440);
  mu[t] <-     mu0 * exp((-E_mu/k)  * (1/temp[t]-1/T0_mu))* (1+ exp(E_muL/k * (1/temp[t]-1/T_muL))  +  exp(E_muH/k * (-1/temp[t]+1/T_muH)));
  beta[t] <- beta0*exp((-E_beta/k) * ((1/temp[t])-(1/T0_beta)))*((1+exp((E_betaH/k)*((-1/temp[t])+(1/T_betaH))))^(-1))
  sigma[t] <-  sigma0*(mass^is)*exp((-E_sigma/k) *  (1/temp[t]-1/T0_sigma)) * (1+  exp(E_sigmaL/k  * (1/temp[t]-1/T_sigmaL)) +  exp(E_sigmaH/k * (-1/temp[t]+1/T_sigmaH)))^(-1); 
  equil.abundance[t] <-  omega0  * exp((-E_omega/k) *  (1/temp[t]-1/T0_omega)) * (1+  exp(E_omegaL/k  * (1/temp[t]-1/T_omegaL)) +  exp(E_omegaH/k * (-1/temp[t]+1/T_omegaH)))^(-1)
  omega[t] <- (equil.abundance[t]*abundance.prop)*spores.per.cluster
  alpha[t] <- (equil.abundance[t]*abundance.prop)*5.12E-6;
  lambda[t] <- (equil.abundance[t]*abundance.prop)*(spores.per.cluster*burst.size)/burst.time
}


#probability of infection 
sigma.prob <- 1 - exp(-sigma*gut.res.time)


#### Parameters ###

### Many parameters will be vectors of length 8 (for each of the 8 temperatures from 10C to 13.5C) ###


# Input from stocks #
# Breakdown of loads from infected stocks is
# Uninfected: 0.5350318
# Infected : 0.4649682

# Input S would be 9 + the probability above * 3
# Input for infecteds would be 3*probabilities above
# Convert to daily by dividing those by 3, aka 9 becomes 3 and the other 3s are removed
inputS = 3 + 0.5350318
inputI = 0.4649682

# Max recruitment rate #
# Can be temperature dependent but probably doesn't need to be if sufficiently high #
recruit = 1.33
  
# Carrying capacity #
# Doesn't need to be temperature dependent #
# This is the mean abundance across all populations throughout the experiment #
kappa = 170

# The next parameters are determined from the MTE models above #
# Index 1001 = 10.0C, Index 1351 = 13.5C #

# Contact rate #
# Temperature dependent, calculate from AmNat paper #
chi.t = c(chi[1001],chi[1051],chi[1101],chi[1151],chi[1201],chi[1251],chi[1301],chi[1351])
  
# Probability of infection#
# Temperature dependent, calculate from AmNat paper #
sigma.t = c(sigma.prob[1001],sigma.prob[1051],sigma.prob[1101],sigma.prob[1151],sigma.prob[1201],sigma.prob[1251],sigma.prob[1301],sigma.prob[1351])
  
# Natural mortality rate #
# Temperature dependent, calculate from PlosBiol paper #

# Temp[101] = 10C, so want 101, 106, 111, 116, 121, 126, 131, 136 to get to 13.5C #

for(i in c(1001,1051,1101,1151,1201,1251,1301,1351)){
  
  mu.sim <- mu[i]
  beta.sim <- beta[i]
  
  U<-1 ; #Starting conditions
  N0<-c(U)
  
  TT<-seq(0.1,400,0.1)   #Set amount of time steps to run in simulation
  step <- 0.1
  
  parms<-c(beta.sim,mu.sim)
  
  survival<-function(t,y,p){
    beta<-p[1]; mu<-p[2];
    U<-y[1];
    
    dU <- -beta*mu^(beta)*(t^(beta-1))*U
    
    list(c(dU))
  }
  
  #run lsoda	
  predictions <-lsoda(N0,TT,survival,parms)	
  
  nam <- paste("predictions",i, sep="")
  assign(nam, predictions)
  
}




#### GET POINT AT WHICH SURVIVAL OF UNINFECTEDS IS 50% ####

surv.10 <- predictions1001[,2]
surv.105 <- predictions1051[,2]
surv.11 <- predictions1101[,2]
surv.115 <- predictions1151[,2]
surv.12 <- predictions1201[,2]
surv.125 <- predictions1251[,2]
surv.13 <- predictions1301[,2]
surv.135 <- predictions1351[,2]

#50% survival between index 1226 and index 1227
# Crosses 1% survival between index 3764 and 3765 #

# 1 / lifespan gives mortality rate #
# Will take lifespan as the lifespan of the average individual #
# Which will be the first timepoint where survival has crossed or equals 50% #

hazard.10 <- 1/predictions1001[which(surv.10 <= 0.5)[[1]],1]
hazard.105 <- 1/predictions1051[which(surv.105 <= 0.5)[[1]],1]
hazard.11 <- 1/predictions1101[which(surv.11 <= 0.5)[[1]],1]
hazard.115 <- 1/predictions1151[which(surv.115 <= 0.5)[[1]],1]
hazard.12 <- 1/predictions1201[which(surv.12 <= 0.5)[[1]],1]
hazard.125 <- 1/predictions1251[which(surv.125 <= 0.5)[[1]],1]
hazard.13 <- 1/predictions1301[which(surv.13 <= 0.5)[[1]],1]
hazard.135 <- 1/predictions1351[which(surv.135 <= 0.5)[[1]],1]

# This is the "new" mu, which is the mortality rate assuming exponential hazard #
mu.t <- as.numeric(c(hazard.10,hazard.105,hazard.11,hazard.115,hazard.12,hazard.125,
             hazard.13,hazard.135))

# Omega #
# Parasites released at death #
omega.t = c(omega[1001],omega[1051],omega[1101],omega[1151],omega[1201],omega[1251],omega[1301],omega[1351])

# Lambda #
# Parasite shedding rate #
lambda.t = c(lambda[1001],lambda[1051],lambda[1101],lambda[1151],lambda[1201],lambda[1251],lambda[1301],lambda[1351])

# Alpha #
# Parasite virulence rate #
alpha.t = c(alpha[1001],alpha[1051],alpha[1101],alpha[1151],alpha[1201],alpha[1251],alpha[1301],alpha[1351])

# Harvesting rate #
# 12 every three days, so 4 per day #
# But per capita would be 0.025 because 4/170 = 0.0235 #
h = 0.0235

# Degradation of dead infecteds rate #
# Just set this to a constant for now #
theta = 0.1
  
# We take 3L out of the 35L every 3 days. Equivalent of 1L per day, or 0.02857143 mortality for spores #
# For now we will assume this is their natural mortality rate (i.e. it swamps out any other environmental mortality)
gamma <-  0.02857143


# organize params here #

parms10 <- c(recruit=recruit,chi.t=chi.t[1],sigma.t=sigma.t[1],mu.t=mu.t[1],alpha.t=alpha.t[1],theta=theta,lambda.t=lambda.t[1],omega.t=omega.t[1],gamma=gamma,h=h,inputS=inputS,inputI=inputI)
parms105 <- c(recruit=recruit,chi.t=chi.t[2],sigma.t=sigma.t[2],mu.t=mu.t[2],alpha.t=alpha.t[2],theta=theta,lambda.t=lambda.t[2],omega.t=omega.t[2],gamma=gamma,h=h,inputS=inputS,inputI=inputI)
parms11 <- c(recruit=recruit,chi.t=chi.t[3],sigma.t=sigma.t[3],mu.t=mu.t[3],alpha.t=alpha.t[3],theta=theta,lambda.t=lambda.t[3],omega.t=omega.t[3],gamma=gamma,h=h,inputS=inputS,inputI=inputI)
parms115 <- c(recruit=recruit,chi.t=chi.t[4],sigma.t=sigma.t[4],mu.t=mu.t[4],alpha.t=alpha.t[4],theta=theta,lambda.t=lambda.t[4],omega.t=omega.t[4],gamma=gamma,h=h,inputS=inputS,inputI=inputI)
parms12 <- c(recruit=recruit,chi.t=chi.t[5],sigma.t=sigma.t[5],mu.t=mu.t[5],alpha.t=alpha.t[5],theta=theta,lambda.t=lambda.t[5],omega.t=omega.t[5],gamma=gamma,h=h,inputS=inputS,inputI=inputI)
parms125 <- c(recruit=recruit,chi.t=chi.t[6],sigma.t=sigma.t[6],mu.t=mu.t[6],alpha.t=alpha.t[6],theta=theta,lambda.t=lambda.t[6],omega.t=omega.t[6],gamma=gamma,h=h,inputS=inputS,inputI=inputI)
parms13 <- c(recruit=recruit,chi.t=chi.t[7],sigma.t=sigma.t[7],mu.t=mu.t[7],alpha.t=alpha.t[7],theta=theta,lambda.t=lambda.t[7],omega.t=omega.t[7],gamma=gamma,h=h,inputS=inputS,inputI=inputI)
parms135 <- c(recruit=recruit,chi.t=chi.t[8],sigma.t=sigma.t[8],mu.t=mu.t[8],alpha.t=alpha.t[8],theta=theta,lambda.t=lambda.t[8],omega.t=omega.t[8],gamma=gamma,h=h,inputS=inputS,inputI=inputI)




# nu is the state=change matrix #

# WOULD BE BETTER IF CAN MAKE THE SAMPLING AND INPUT PARAMS DETERMINISTIC

nu <- matrix(c(+1,-1,-1,0,0,0,0,0,0,-1,0,+1,0,
               0,+1,0,-1,-1,0,0,0,0,0,-1,0,+1,
               0,0,0,+1,+1,-1,0,0,0,0,0,0,0,
               0,0,0,0,0,0,+1,+1,-1,0,0,0,0),
             nrow=4,byrow=TRUE)


# a is the propensity functions which need to be derived from the model #
# these are the 14 reactions in the system #

a <- c("recruit",   # 1 - density-independent recruitment
       "chi.t*sigma.t*S*E",                # 3 - transmission
       "mu.t*S",                     # 4 - natural death of susceptible
       "mu.t*I",                     # 5 - natural death of infected
       "alpha.t*I",                  # 6 - parasite-induced death of infected
       "theta*D",                  # 7 - decay and loss of infected corpse
       "lambda.t*I",                 # 8 -parasite spore shed by infected
       "omega.t*theta*D",            # 9 - parasite spore shed by corpse
       "gamma*E",                  # 10 - parasite death
       "h*(S+I)*(S/(S+I))",        # 11 - sampling uninfected
       "h*(S+I)*(I/(S+I))",        # 12 - sampling infected
       "inputS",                  # 13 - inputting susceptible
       "inputI")                  # 14 - inputting infected



out10 <- NULL
out105 <- NULL
out11 <- NULL
out115 <- NULL
out12 <- NULL
out125 <- NULL
out13 <- NULL
out135 <- NULL


x10 <- NULL
x105 <- NULL
x11 <- NULL
x115 <- NULL
x12 <- NULL
x125 <- NULL
x13 <- NULL
x135 <- NULL



### sim ###


number.of.sims = 4

# Initial conditions
x10 <- c(S=169,I=0,D=0,E=0)

set.seed <- (1234)


for(i in 1:number.of.sims){
out10[[i]] <- ssa(x10, a, nu, parms10, tf = 15, method = "ETL",
                   simName = "SIDE Model",
                   verbose = TRUE,
                   consoleInterval = Inf,
                   tau=0.001)
}


# 10.5 degrees #
for(i in 1:number.of.sims){
  x105[[i]] <- c(S=out10[[i]]$data[15003,2],I=out10[[i]]$data[15003,3],D=out10[[i]]$data[15003,4],E=out10[[i]]$data[15003,5])
}

for(i in 1:number.of.sims){
  out105[[i]]<- ssa(x105[[i]], a, nu, parms105, tf = 15, method = "ETL",
                    simName = "SIDE Model",
                    verbose = TRUE,
                    consoleInterval = Inf,
                    tau=0.001)
}


# 11 degrees #
for(i in 1:number.of.sims){
  x11[[i]] <- c(S=out105[[i]]$data[15003,2],I=out105[[i]]$data[15003,3],D=out105[[i]]$data[15003,4],E=out105[[i]]$data[15003,5])
}


for(i in 1:number.of.sims){
  out11[[i]]<- ssa(x11[[i]], a, nu, parms11, tf = 15, method = "ETL",
                   simName = "SIDE Model",
                   verbose = TRUE,
                   consoleInterval = Inf,
                   tau=0.001)
}

# 11.5 degrees #
for(i in 1:number.of.sims){
  x115[[i]] <- c(S=out11[[i]]$data[15003,2],I=out11[[i]]$data[15003,3],D=out11[[i]]$data[15003,4],E=out11[[i]]$data[15003,5])
}

for(i in 1:number.of.sims){
  out115[[i]]<- ssa(x115[[i]], a, nu, parms115, tf = 15, method = "ETL",
                    simName = "SIDE Model",
                    verbose = TRUE,
                    consoleInterval = Inf,
                    tau=0.001)
}


# 12 degrees #
for(i in 1:number.of.sims){
  x12[[i]] <- c(S=out115[[i]]$data[15003,2],I=out115[[i]]$data[15003,3],D=out115[[i]]$data[15003,4],E=out115[[i]]$data[15003,5])
}

for(i in 1:number.of.sims){
  out12[[i]]<- ssa(x12[[i]], a, nu, parms12, tf = 15, method = "ETL",
                   simName = "SIDE Model",
                   verbose = TRUE,
                   consoleInterval = Inf,
                   tau=0.001)
}



# 12.5 degrees #
for(i in 1:number.of.sims){
  x125[[i]] <- c(S=out12[[i]]$data[15003,2],I=out12[[i]]$data[15003,3],D=out12[[i]]$data[15003,4],E=out12[[i]]$data[15003,5])
}

for(i in 1:number.of.sims){
  out125[[i]]<- ssa(x125[[i]], a, nu, parms125, tf = 15, method = "ETL",
                    simName = "SIDE Model",
                    verbose = TRUE,
                    consoleInterval = Inf,
                    tau=0.001)
}


# 13 degrees #
for(i in 1:number.of.sims){
  x13[[i]] <- c(S=out125[[i]]$data[15003,2],I=out125[[i]]$data[15003,3],D=out125[[i]]$data[15003,4],E=out125[[i]]$data[15003,5])
}

for(i in 1:number.of.sims){
  out13[[i]]<- ssa(x13[[i]], a, nu, parms13, tf = 15, method = "ETL",
                   simName = "SIDE Model",
                   verbose = TRUE,
                   consoleInterval = Inf,
                   tau=0.001)
}

# 13.5 degrees #
for(i in 1:number.of.sims){
  x135[[i]] <- c(S=out13[[i]]$data[15003,2],I=out13[[i]]$data[15003,3],D=out13[[i]]$data[15003,4],E=out13[[i]]$data[15003,5])
}

for(i in 1:number.of.sims){
  out135[[i]]<- ssa(x135[[i]], a, nu, parms135, tf = 15, method = "ETL",
                    simName = "SIDE Model",
                    verbose = TRUE,
                    consoleInterval = Inf,
                    tau=0.001)
}






# Organize it all together #
warming <- NULL

for(i in 1:number.of.sims){
  warming[[i]] <- rbind(out10[[i]]$data,out105[[i]]$data,out11[[i]]$data,out115[[i]]$data,
                         out12[[i]]$data,out125[[i]]$data,out13[[i]]$data,
                         out135[[i]]$data)
  
  warming[[i]][,1] <- seq(0.000,120.023,0.001)
  colnames(warming[[i]]) <- c("Time","S","I","D","E")
}



#### Control populations that stay at 10 degrees ###


constant <- NULL

# 10 degrees constant #
for(i in 1:number.of.sims){
  constant[[i]] <- ssa(x10, a, nu, parms10, tf = 120, method = "ETL",
                       simName = "SIDE Model",
                       verbose = TRUE,
                       consoleInterval = Inf,
                       tau=0.001)
}





#### constant and warming on same graph ####

quartz()
par(mar=c(5,5,3,2))
plot(constant[[1]]$data[,1],(constant[[1]]$data[,3]/(constant[[1]]$data[,2]+constant[[1]]$data[,3])),type='n',xlim=c(0,120),ylim=c(0,1),xlab="Time (days)",ylab="Prevalence",
     main=NA,cex.lab=2,cex.main=2)

for(i in 1:length(constant)){
  lines(constant[[i]]$data[,1],(constant[[i]]$data[,3]/(constant[[i]]$data[,2]+constant[[i]]$data[,3])),col="deepskyblue",lwd=2)
}

for(i in 1:length(warming)){
  lines(warming[[i]][,1],(warming[[i]][,3]/(warming[[i]][,2]+warming[[i]][,3])),col="black",lwd=2)
}







# LOAD IMAGE #


# Get 95% quantiles of prevalence to plot with epidemics data #


warming.prevalences <- matrix(nrow=120024,ncol=4)
constant.prevalences <- matrix(nrow=120002,ncol=4)


for(i in 1:4){
warming.prevalences[,i] <- warming[[i]][,3]/(warming[[i]][,3]+warming[[i]][,2])
}
for(i in 1:4){
constant.prevalences[,i] <- as.numeric(constant[[i]]$data[,3])/(as.numeric(constant[[i]]$data[,3])+as.numeric(constant[[i]]$data[,2]))
}


# Get confidence intervals and medians for MTE predictions #
upper.warming.prev <- c()
lower.warming.prev <- c()
median.warming.prev <- c()

upper.constant.prev <- c()
lower.constant.prev <- c()
median.constant.prev <- c()

for(i in 1:dim(warming.prevalences)[[1]]){
  upper.warming.prev[i] <- quantile(warming.prevalences[i,], probs = 0.975)
}

for(i in 1:dim(warming.prevalences)[[1]]){
  lower.warming.prev[i] <- quantile(warming.prevalences[i,], probs = 0.025)
}

for(i in 1:dim(warming.prevalences)[[1]]){
  median.warming.prev[i] <- median(warming.prevalences[i,])
}



for(i in 1:dim(constant.prevalences)[[1]]){
  upper.constant.prev[i] <- quantile(constant.prevalences[i,], probs = 0.975)
}

for(i in 1:dim(constant.prevalences)[[1]]){
  lower.constant.prev[i] <- quantile(constant.prevalences[i,], probs = 0.025)
}

for(i in 1:dim(constant.prevalences)[[1]]){
  median.constant.prev[i] <- median(constant.prevalences[i,])
}





require(gdata)
require(doBy)
require(ecp)
require(strucchange)


data <- read.xls("~/Desktop/Epidemics Project/epidemics_data.xlsx", header=TRUE)

data.sum <- summaryBy(status ~ population+treatment+day, data = data, FUN = sum)

sub1 <- subset(data.sum, population == 1)
sub2 <- subset(data.sum, population == 2)
sub3 <- subset(data.sum, population == 3)
sub4 <- subset(data.sum, population == 4)
sub5 <- subset(data.sum, population == 5)
sub6 <- subset(data.sum, population == 6)
sub7 <- subset(data.sum, population == 7)
sub8 <- subset(data.sum, population == 8)

warming.time.axis <- c(seq(0.001,120.024,0.001))
constant.time.axis <- c(seq(0.001,120.002,0.001))


###### Plot all 8 separately ########

quartz()
par(mfrow=c(2,4))
par(mar=c(1,1,1,1))
par(oma=c(6,6,0,0))
plot(sub1$status.sum/12 ~ sub1$day, ylim=c(0,1), xlim=c(0,120), col="blue", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = max(sub1$day)/2, y = 0.9, labels = "10C", cex=1.5)
lines(constant.time.axis, lower.constant.prev,col="black",lwd=2)
lines(constant.time.axis, upper.constant.prev,col="black",lwd=2)
lines(constant.time.axis, median.constant.prev, col="black", lwd=3)

plot(sub2$status.sum/12 ~ sub2$day, ylim=c(0,1), xlim=c(0,120),col="blue", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = max(sub1$day)/2, y = 0.9, labels = "10C", cex=1.5)
lines(constant.time.axis, lower.constant.prev,col="black",lwd=2)
lines(constant.time.axis, upper.constant.prev,col="black",lwd=2)
lines(constant.time.axis, median.constant.prev, col="black", lwd=3)

plot(sub3$status.sum/12 ~ sub3$day, ylim=c(0,1), xlim=c(0,120),col="blue", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = max(sub1$day)/2, y = 0.9, labels = "10C", cex=1.5)
lines(constant.time.axis, lower.constant.prev,col="black",lwd=2)
lines(constant.time.axis, upper.constant.prev,col="black",lwd=2)
lines(constant.time.axis, median.constant.prev, col="black", lwd=3)

plot(sub4$status.sum/12 ~ sub4$day, ylim=c(0,1), xlim=c(0,120),col="blue", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = max(sub1$day)/2, y = 0.9, labels = "10C", cex=1.5)
lines(constant.time.axis, lower.constant.prev,col="black",lwd=2)
lines(constant.time.axis, upper.constant.prev,col="black",lwd=2)
lines(constant.time.axis, median.constant.prev, col="black", lwd=3)

plot(sub5$status.sum/12 ~ sub5$day, ylim=c(0,1), xlim=c(0,120),col="red", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = c(7.5,22.5,37.5,52.5,67.5,82.5,97.5,112.5),y = 0.9, labels = c("10","10.5","11","11.5","12","12.5","13","13.5"), cex=1.5)
lines(warming.time.axis, lower.warming.prev,col="black",lwd=2)
lines(warming.time.axis, upper.warming.prev,col="black",lwd=2)
lines(warming.time.axis, median.warming.prev, col="black", lwd=3)
abline(v=c(15,30,45,60,75,90,105), lty=3)


plot(sub6$status.sum/12 ~ sub6$day, ylim=c(0,1), xlim=c(0,120),col="red", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = c(7.5,22.5,37.5,52.5,67.5,82.5,97.5,112.5),y = 0.9, labels = c("10","10.5","11","11.5","12","12.5","13","13.5"), cex=1.5)
lines(warming.time.axis, lower.warming.prev,col="black",lwd=2)
lines(warming.time.axis, upper.warming.prev,col="black",lwd=2)
lines(warming.time.axis, median.warming.prev, col="black", lwd=3)
abline(v=c(15,30,45,60,75,90,105), lty=3)



plot(sub7$status.sum/12 ~ sub7$day, ylim=c(0,1), xlim=c(0,120),col="red", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = c(7.5,22.5,37.5,52.5,67.5,82.5,97.5,112.5),y = 0.9, labels = c("10","10.5","11","11.5","12","12.5","13","13.5"), cex=1.5)
lines(warming.time.axis, lower.warming.prev,col="black",lwd=2)
lines(warming.time.axis, upper.warming.prev,col="black",lwd=2)
lines(warming.time.axis, median.warming.prev, col="black", lwd=3)
abline(v=c(15,30,45,60,75,90,105), lty=3)

plot(sub8$status.sum/12 ~ sub8$day, ylim=c(0,1), xlim=c(0,120),col="red", type="l", lty=1,
     xlab=NA, ylab=NA, lwd=4, cex.lab=1.5, pch=19)
text(x = c(7.5,22.5,37.5,52.5,67.5,82.5,97.5,112.5),y = 0.9, labels = c("10","10.5","11","11.5","12","12.5","13","13.5"), cex=1.5)
lines(warming.time.axis, lower.warming.prev,col="black",lwd=2)
lines(warming.time.axis, upper.warming.prev,col="black",lwd=2)
lines(warming.time.axis, median.warming.prev, col="black", lwd=3)
abline(v=c(15,30,45,60,75,90,105), lty=3)


mtext("Prevalence",side=2, cex=3, line = 2.5, outer=TRUE)
mtext("Time (days)",side=1, cex=3, line = 3.5, outer=TRUE)





