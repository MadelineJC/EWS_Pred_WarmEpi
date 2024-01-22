require(scales)
require(evobiR)
require(gdata)
require(doBy)
require(ecp)
require(strucchange)
require(TTR)
require(moments)

# change R0 ablines (check other code first for temperature again)
# change legend coordinates
# calculate 2.5 and 97.5% quantiles
# plot CIs

# Model predicts R0 > 1.0 at 11.71 degrees #
# So for warming sims this is first when the temp is 12.0C #
# And 12.0C starts at the halfway point of the experiment #


num.of.sims <- 1000

constant1 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.1.csv")
constant2 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.2.csv")
constant3 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.3.csv")
constant4 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.4.csv")
constant5 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.5.csv")
constant6 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.6.csv")
constant7 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.7.csv")
constant8 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.8.csv")
constant9 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.9.csv")
constant10 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/constant.incidence.last.120.10.csv")

constant.incidence.daily <- cbind(constant1[,2:101],constant2[,2:101],constant3[,2:101],constant4[,2:101],
                                  constant5[,2:101],constant6[,2:101],constant7[,2:101],
                                  constant8[,2:101],constant9[,2:101],constant10[,2:101])

warning1 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.1.csv")
warning2 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.2.csv")
warning3 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.3.csv")
warning4 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.4.csv")
warning5 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.5.csv")
warning6 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.6.csv")
warning7 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.7.csv")
warning8 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.8.csv")
warning9 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.9.csv")
warning10 <- read.csv("Kirk_Data/Sim_Pops_EWS_Analyses/warning.incidence.last.120.10.csv")


warming.incidence.daily <- cbind(warming1[,2:101],warming2[,2:101],warming3[,2:101],warming4[,2:101],
                                 warming5[,2:101],warming6[,2:101],warming7[,2:101],
                                 warming8[,2:101],warming9[,2:101],warming10[,2:101])



# PLOTTING 1000 SIMS OF BOTH TYPES IS A LOT FOR NOW, SO PLOT 1:100 OF EACH FOR NOW #
quartz()
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
plot(constant.incidence.daily[,1], col="blue", type="n",ylim=c(0,100),xlim=c(0,120),
     xlab="Time",ylab="Number of Infecteds",cex.lab=1.5)
for(i in 1:100){
  lines(constant.incidence.daily[,i], col=alpha("blue",0.2), lwd=3)
  lines(warming.incidence.daily[,i], col=alpha("red",0.2), lwd=3)
}
abline(v=60, lty=2, col="black", lwd=2)
legend(0,100,legend=c("Warming",expression('Constant 10'~degree~'C'),expression('R'['0']*' crosses 1')),lwd=c(3,3,2),col=c("red","blue","black"),lty=c(1,1,2))


## look at incidence right before R0 crosses ##
quartz()
par(mfrow=c(1,1))
par(mar=c(5,5,1,1))
plot(constant.incidence.daily[,1], col="blue", type="n",ylim=c(0,40),xlim=c(40,70),
     xlab="Time",ylab="Number of Infecteds",cex.lab=1.5)
for(i in 1:100){
  lines(constant.incidence.daily[,i], col=alpha("blue",0.2), lwd=3)
  lines(warming.incidence.daily[,i], col=alpha("red",0.2), lwd=3)
}
abline(v=60, lty=2, col="black", lwd=2)
legend(0,100,legend=c("Warming",expression('Constant 10'~degree~'C'),expression('R'['0']*' crosses 1')),lwd=c(3,3,2),col=c("red","blue","black"),lty=c(1,1,2))








### WINDOW SIZE OF 5 DAYS ####

# FOR LOOPS BELOW COMPUTE SD, MEAN, SKEWNESS, AND KURTOSIS #

# Sliding window function doesn't calculate the window to include the very last index, but this is okay since this is long after the epidemic transition should have occurred #

warming.sims.sd.5 <- matrix(ncol=num.of.sims,nrow=115)
warming.sims.mean.5 <- matrix(ncol=num.of.sims,nrow=115)
warming.sims.skewness.5 <- matrix(ncol=num.of.sims,nrow=115)
warming.sims.kurtosis.5 <- matrix(ncol=num.of.sims,nrow=115)


for(i in 1:num.of.sims){
  warming.sims.sd.5[,i] <- SlidingWindow(FUN="sd",data=warming.incidence.daily[,i], window=5,step=1)
  warming.sims.mean.5[,i] <- SlidingWindow(FUN="mean",data=warming.incidence.daily[,i], window=5,step=1)
  warming.sims.skewness.5[,i] <- SlidingWindow(FUN="skewness",data=warming.incidence.daily[,i], window=5,step=1)
  warming.sims.kurtosis.5[,i] <- SlidingWindow(FUN="kurtosis",data=warming.incidence.daily[,i], window=5,step=1)
  
}


constant.sims.sd.5 <- matrix(ncol=num.of.sims,nrow=115)
constant.sims.mean.5 <- matrix(ncol=num.of.sims,nrow=115)
constant.sims.skewness.5 <- matrix(ncol=num.of.sims,nrow=115)
constant.sims.kurtosis.5 <- matrix(ncol=num.of.sims,nrow=115)

for(i in 1:num.of.sims){
  constant.sims.sd.5[,i] <- SlidingWindow(FUN="sd",data=constant.incidence.daily[,i], window=5,step=1)
  constant.sims.mean.5[,i] <- SlidingWindow(FUN="mean",data=constant.incidence.daily[,i], window=5,step=1)
  constant.sims.skewness.5[,i] <- SlidingWindow(FUN="skewness",data=constant.incidence.daily[,i], window=5,step=1)
  constant.sims.kurtosis.5[,i] <- SlidingWindow(FUN="kurtosis",data=constant.incidence.daily[,i], window=5,step=1)
}



## VARIANCE ##
# sd * sd # # note that using the square (^) function multiples the wrong parts of the matrix together #
warming.sims.variance.5 <- warming.sims.sd.5*warming.sims.sd.5
constant.sims.variance.5 <- constant.sims.sd.5*constant.sims.sd.5


## COEFFICIENT OF VARIATION ##
# coefficient of variation = sd / mean #
warming.sims.cv.5 <- warming.sims.sd.5/warming.sims.mean.5
constant.sims.cv.5 <- constant.sims.sd.5/constant.sims.mean.5



# INDEX OF DISPERSION #
# index of dispersion = variance / mean #
warming.sims.id.5 <- warming.sims.variance.5/warming.sims.mean.5
constant.sims.id.5 <- constant.sims.variance.5/constant.sims.mean.5


# FIRST DIFFERENCE VARIANCE #
# Variance(t) - Variance(t-1)
# Output will be 1 index shorter than the rest of the metrics, as the final index will be NA #

warming.sims.first.diff.variance.5 <- matrix(ncol=num.of.sims,nrow=115)
constant.sims.first.diff.variance.5 <- matrix(ncol=num.of.sims,nrow=115)

for(i in 1:num.of.sims){
  for(j in 2:115){
    warming.sims.first.diff.variance.5[j-1,i] <- warming.sims.variance.5[j,i] - warming.sims.variance.5[j-1,i]
    constant.sims.first.diff.variance.5[j-1,i] <- constant.sims.variance.5[j,i] - constant.sims.variance.5[j-1,i]
  }
}


### RUN THIS FOR BOTH AUTOCORRELATION AND AUTOCOVARIANCE ###
# using acf(segment)[[1]][2] pulls out the lag-1 acf

# Since sliding window doesn't calculate metrics to include the very last index, we won't here for cov or cor either #

constant.sims.cor.5 <- matrix(ncol=num.of.sims,nrow=120)
constant.sims.cov.5 <- matrix(ncol=num.of.sims,nrow=120)

for(j in 1:num.of.sims){
  for(i in 1:115){
    window.start <- i
    window.end <- i+4
    
    segment <- constant.incidence.daily[c(window.start:window.end),j]
    
    constant.sims.cor.5[i,j] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
    constant.sims.cov.5[i,j] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
    
  }
}



warming.sims.cor.5 <- matrix(ncol=num.of.sims,nrow=120)
warming.sims.cov.5 <- matrix(ncol=num.of.sims,nrow=120)

for(j in 1:num.of.sims){
  for(i in 1:115){
    window.start <- i
    window.end <- i+4
    
    segment <- warming.incidence.daily[c(window.start:window.end),j]
    
    warming.sims.cor.5[i,j] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
    warming.sims.cov.5[i,j] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
    
  }
}

#### DECAY TIME ###
# -time step / ln(autocorrelation) #
# This will create some NaNs, so use supress warnings function to surpress them #


warming.sims.decay.time.5 <- matrix(ncol=num.of.sims,nrow=55)
constant.sims.decay.time.5 <- matrix(ncol=num.of.sims,nrow=55)

for(i in 1:num.of.sims){
  for(j in 1:55){
    suppressWarnings(warming.sims.decay.time.5[j,i] <- -j/log(warming.sims.cor.5[j,i]))
    suppressWarnings(constant.sims.decay.time.5[j,i] <- -j/log(constant.sims.cor.5[j,i]))
  }
}

################################################################################################################################
################################################################################################################################
################################################################################################################################







### WINDOW SIZE OF 15 DAYS ####

# for LOOPS BELOW COMPUTE SD, MEAN, SKEWNESS, AND KURTOSIS #

warming.sims.sd.15 <- matrix(ncol=num.of.sims,nrow=105)
warming.sims.mean.15 <- matrix(ncol=num.of.sims,nrow=105)
warming.sims.skewness.15 <- matrix(ncol=num.of.sims,nrow=105)
warming.sims.kurtosis.15 <- matrix(ncol=num.of.sims,nrow=105)


for(i in 1:num.of.sims){
  warming.sims.sd.15[,i] <- SlidingWindow(FUN="sd",data=warming.incidence.daily[,i], window=15,step=1)
  warming.sims.mean.15[,i] <- SlidingWindow(FUN="mean",data=warming.incidence.daily[,i], window=15,step=1)
  warming.sims.skewness.15[,i] <- SlidingWindow(FUN="skewness",data=warming.incidence.daily[,i], window=15,step=1)
  warming.sims.kurtosis.15[,i] <- SlidingWindow(FUN="kurtosis",data=warming.incidence.daily[,i], window=15,step=1)
  
}


constant.sims.sd.15 <- matrix(ncol=num.of.sims,nrow=105)
constant.sims.mean.15 <- matrix(ncol=num.of.sims,nrow=105)
constant.sims.skewness.15 <- matrix(ncol=num.of.sims,nrow=105)
constant.sims.kurtosis.15 <- matrix(ncol=num.of.sims,nrow=105)

for(i in 1:num.of.sims){
  constant.sims.sd.15[,i] <- SlidingWindow(FUN="sd",data=constant.incidence.daily[,i], window=15,step=1)
  constant.sims.mean.15[,i] <- SlidingWindow(FUN="mean",data=constant.incidence.daily[,i], window=15,step=1)
  constant.sims.skewness.15[,i] <- SlidingWindow(FUN="skewness",data=constant.incidence.daily[,i], window=15,step=1)
  constant.sims.kurtosis.15[,i] <- SlidingWindow(FUN="kurtosis",data=constant.incidence.daily[,i], window=15,step=1)
}


# COEFFICIENT OF VARIATION #
warming.sims.cv.15 <- warming.sims.sd.15/warming.sims.mean.15
constant.sims.cv.15 <- constant.sims.sd.15/constant.sims.mean.15


## VARIANCE ##
warming.sims.variance.15 <- warming.sims.sd.15*warming.sims.sd.15
constant.sims.variance.15 <- constant.sims.sd.15*constant.sims.sd.15


# INDEX OF DISPERSION #
# index of dispersion = variance / mean #
warming.sims.id.15 <- warming.sims.variance.15/warming.sims.mean.15
constant.sims.id.15 <- constant.sims.variance.15/constant.sims.mean.15


# FIRST DIFFERENCE VARIANCE #
# Variance(t) - Variance(t-1)
# Output will be 1 index shorter than the rest of the metrics, as the final index will be NA #

warming.sims.first.diff.variance.15 <- matrix(ncol=num.of.sims,nrow=105)
constant.sims.first.diff.variance.15 <- matrix(ncol=num.of.sims,nrow=105)

for(i in 1:num.of.sims){
  for(j in 2:105){
    warming.sims.first.diff.variance.15[j-1,i] <- warming.sims.variance.15[j,i] - warming.sims.variance.15[j-1,i]
    constant.sims.first.diff.variance.15[j-1,i] <- constant.sims.variance.15[j,i] - constant.sims.variance.15[j-1,i]
  }
}




# using acf(segment)[[1]][2] pulls out the lag-1 acf

# Since sliding window doesn't calculate metrics to include the very last index, we won't here for cov or cor either #

constant.sims.cor.15 <- matrix(ncol=num.of.sims,nrow=120)
constant.sims.cov.15 <- matrix(ncol=num.of.sims,nrow=120)


for(j in 1:num.of.sims){
  for(i in 1:105){
    window.start <- i
    window.end <- i+14
    
    segment <- constant.incidence.daily[c(window.start:window.end),j]
    
    constant.sims.cor.15[i,j] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
    constant.sims.cov.15[i,j] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
    
  }
}


warming.sims.cor.15 <- matrix(ncol=num.of.sims,nrow=120)
warming.sims.cov.15 <- matrix(ncol=num.of.sims,nrow=120)


for(j in 1:num.of.sims){
  for(i in 1:105){
    window.start <- i
    window.end <- i+14
    
    segment <- warming.incidence.daily[c(window.start:window.end),j]
    
    warming.sims.cor.15[i,j] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
    warming.sims.cov.15[i,j] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
    
  }
}


#### DECAY TIME ###
# -time step / ln(autocorrelation) #
# This will create some NaNs, so use supress warnings function to surpress them #

warming.sims.decay.time.15 <- matrix(ncol=num.of.sims,nrow=45)
constant.sims.decay.time.15 <- matrix(ncol=num.of.sims,nrow=45)

for(i in 1:num.of.sims){
  for(j in 1:45){
    suppressWarnings(warming.sims.decay.time.15[j,i] <- -j/log(warming.sims.cor.15[j,i]))
    suppressWarnings(constant.sims.decay.time.15[j,i] <- -j/log(constant.sims.cor.15[j,i]))
  }
}





################################################################################################################################
################################################################################################################################
################################################################################################################################




#### Now repeat everything with window size of 30 days #

warming.sims.sd.30 <- matrix(ncol=num.of.sims,nrow=90)
warming.sims.mean.30 <- matrix(ncol=num.of.sims,nrow=90)
warming.sims.skewness.30 <- matrix(ncol=num.of.sims,nrow=90)
warming.sims.kurtosis.30 <- matrix(ncol=num.of.sims,nrow=90)


for(i in 1:num.of.sims){
  warming.sims.sd.30[,i] <- SlidingWindow(FUN="sd",data=warming.incidence.daily[,i], window=30,step=1)
  warming.sims.mean.30[,i] <- SlidingWindow(FUN="mean",data=warming.incidence.daily[,i], window=30,step=1)
  warming.sims.skewness.30[,i] <- SlidingWindow(FUN="skewness",data=warming.incidence.daily[,i], window=30,step=1)
  warming.sims.kurtosis.30[,i] <- SlidingWindow(FUN="kurtosis",data=warming.incidence.daily[,i], window=30,step=1)
  
}


constant.sims.sd.30 <- matrix(ncol=num.of.sims,nrow=90)
constant.sims.mean.30 <- matrix(ncol=num.of.sims,nrow=90)
constant.sims.skewness.30 <- matrix(ncol=num.of.sims,nrow=90)
constant.sims.kurtosis.30 <- matrix(ncol=num.of.sims,nrow=90)

for(i in 1:num.of.sims){
  constant.sims.sd.30[,i] <- SlidingWindow(FUN="sd",data=constant.incidence.daily[,i], window=30,step=1)
  constant.sims.mean.30[,i] <- SlidingWindow(FUN="mean",data=constant.incidence.daily[,i], window=30,step=1)
  constant.sims.skewness.30[,i] <- SlidingWindow(FUN="skewness",data=constant.incidence.daily[,i], window=30,step=1)
  constant.sims.kurtosis.30[,i] <- SlidingWindow(FUN="kurtosis",data=constant.incidence.daily[,i], window=30,step=1)
}





### COEFFICIENT OF VARIATION ###
warming.sims.cv.30 <- warming.sims.sd.30/warming.sims.mean.30
constant.sims.cv.30 <- constant.sims.sd.30/constant.sims.mean.30


### VARIANCE ###
warming.sims.variance.30 <- warming.sims.sd.30*warming.sims.sd.30
constant.sims.variance.30 <- constant.sims.sd.30*constant.sims.sd.30


# INDEX OF DISPERSION #
# index of dispersion = variance / mean #
warming.sims.id.30 <- warming.sims.variance.30/warming.sims.mean.30
constant.sims.id.30 <- constant.sims.variance.30/constant.sims.mean.30


# FIRST DIFFERENCE VARIANCE #
# Variance(t) - Variance(t-1)
# Output will be 1 index shorter than the rest of the metrics, as the final index will be NA #

warming.sims.first.diff.variance.30 <- matrix(ncol=num.of.sims,nrow=90)
constant.sims.first.diff.variance.30 <- matrix(ncol=num.of.sims,nrow=90)

for(i in 1:num.of.sims){
  for(j in 2:90){
    warming.sims.first.diff.variance.30[j-1,i] <- warming.sims.variance.30[j,i] - warming.sims.variance.30[j-1,i]
    constant.sims.first.diff.variance.30[j-1,i] <- constant.sims.variance.30[j,i] - constant.sims.variance.30[j-1,i]
  }
}



# Since sliding window doesn't calculate metrics to include the very last index, we won't here for cov or cor either #

# using acf(segment)[[1]][2] pulls out the lag-1 acf

constant.sims.cor.30 <- matrix(ncol=num.of.sims,nrow=120)
constant.sims.cov.30 <- matrix(ncol=num.of.sims,nrow=120)

for(j in 1:num.of.sims){
  for(i in 1:90){
    window.start <- i
    window.end <- i+29
    
    segment <- constant.incidence.daily[c(window.start:window.end),j]
    
    constant.sims.cor.30[i,j] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
    constant.sims.cov.30[i,j] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
    
  }
}


warming.sims.cor.30 <- matrix(ncol=num.of.sims,nrow=120)
warming.sims.cov.30 <- matrix(ncol=num.of.sims,nrow=120)


for(j in 1:num.of.sims){
  for(i in 1:90){
    window.start <- i
    window.end <- i+29
    
    segment <- warming.incidence.daily[c(window.start:window.end),j]
    
    warming.sims.cor.30[i,j] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
    warming.sims.cov.30[i,j] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
    
  }
}


#### DECAY TIME ###
# -time step / ln(autocorrelation) #
# This will create some NaNs, so use supress warnings function to surpress them #


warming.sims.decay.time.30 <- matrix(ncol=num.of.sims,nrow=30)
constant.sims.decay.time.30 <- matrix(ncol=num.of.sims,nrow=30)

for(i in 1:num.of.sims){
  for(j in 1:30){
    suppressWarnings(warming.sims.decay.time.30[j,i] <- -j/log(warming.sims.cor.30[j,i]))
    suppressWarnings(constant.sims.decay.time.30[j,i] <- -j/log(constant.sims.cor.30[j,i]))
  }
}


################################################################################################################################
################################################################################################################################
################################################################################################################################








############## TO RUN AUC TESTS (BRETT ET AL. 2018, PLOS COMP. BIOL) TO SEE IF EWS ARE EFFECTIVE, WE NEED TO ONLY LOOK AT THE METRICS BEFORE THE TRANSITIION IS REACHED ######

### Window size of 5 first ###

warming.sims.variance.5 <- warming.sims.variance.5[1:55,]
warming.sims.mean.5 <- warming.sims.mean.5[1:55,]
warming.sims.cv.5 <- warming.sims.cv.5[1:55,]
warming.sims.skewness.5 <-  warming.sims.skewness.5[1:55,]
warming.sims.kurtosis.5 <-  warming.sims.kurtosis.5[1:55,]
warming.sims.id.5 <-  warming.sims.id.5[1:55,]
warming.sims.cor.5 <-  warming.sims.cor.5[1:55,]
warming.sims.cov.5 <-  warming.sims.cov.5[1:55,]
warming.sims.first.diff.variance.5 <-  warming.sims.first.diff.variance.5[1:55,]
warming.sims.decay.time.5 <-  warming.sims.decay.time.5[1:55,]


constant.sims.variance.5 <- constant.sims.variance.5[1:55,]
constant.sims.mean.5 <- constant.sims.mean.5[1:55,]
constant.sims.cv.5 <- constant.sims.cv.5[1:55,]
constant.sims.skewness.5 <-  constant.sims.skewness.5[1:55,]
constant.sims.kurtosis.5 <-  constant.sims.kurtosis.5[1:55,]
constant.sims.id.5 <-  constant.sims.id.5[1:55,]
constant.sims.cor.5 <-  constant.sims.cor.5[1:55,]
constant.sims.cov.5 <-  constant.sims.cov.5[1:55,]
constant.sims.first.diff.variance.5 <-  constant.sims.first.diff.variance.5[1:55,]
constant.sims.decay.time.5 <-  constant.sims.decay.time.5[1:55,]


### Time series for window of 5 ###
time.5 <- c(1:55)

warming.sims.variance.5.tau <- c()
warming.sims.mean.5.tau <- c()
warming.sims.cv.5.tau <- c()
warming.sims.skewness.5.tau <-  c()
warming.sims.kurtosis.5.tau <- c()
warming.sims.id.5.tau <-  c()
warming.sims.cor.5.tau <-  c()
warming.sims.cov.5.tau <-  c()
warming.sims.first.diff.variance.5.tau <-  c()
warming.sims.decay.time.5.tau <-  c()


constant.sims.variance.5.tau <- c()
constant.sims.mean.5.tau <- c()
constant.sims.cv.5.tau <- c()
constant.sims.skewness.5.tau <-  c()
constant.sims.kurtosis.5.tau <- c()
constant.sims.id.5.tau <- c()
constant.sims.cor.5.tau <- c()
constant.sims.cov.5.tau <- c()
constant.sims.first.diff.variance.5.tau <-  c()
constant.sims.decay.time.5.tau <- c()



# Calculate Kendall's tau for each metric vs the time series #
for(i in 1:num.of.sims){
  warming.sims.variance.5.tau[i] <- cor(warming.sims.variance.5[,i], time.5, method="kendall")
  warming.sims.mean.5.tau[i] <- cor(warming.sims.mean.5[,i], time.5, method="kendall")
  warming.sims.cv.5.tau[i] <- cor(warming.sims.cv.5[,i], time.5, method="kendall")
  warming.sims.skewness.5.tau[i] <- cor(warming.sims.skewness.5[,i], time.5, method="kendall")
  warming.sims.kurtosis.5.tau[i] <- cor(warming.sims.kurtosis.5[,i], time.5, method="kendall")
  warming.sims.id.5.tau[i] <- cor(warming.sims.id.5[,i], time.5, method="kendall")
  warming.sims.cor.5.tau[i] <- cor(warming.sims.cor.5[,i], time.5, method="kendall")
  warming.sims.cov.5.tau[i] <- cor(warming.sims.cov.5[,i], time.5, method="kendall")
  warming.sims.first.diff.variance.5.tau[i] <- cor(warming.sims.first.diff.variance.5[,i], time.5, method="kendall")
  warming.sims.decay.time.5.tau[i] <- cor(warming.sims.decay.time.5[,i], time.5, method="kendall")
  
  constant.sims.variance.5.tau[i] <- cor(constant.sims.variance.5[,i], time.5, method="kendall")
  constant.sims.mean.5.tau[i] <- cor(constant.sims.mean.5[,i], time.5, method="kendall")
  constant.sims.cv.5.tau[i] <- cor(constant.sims.cv.5[,i], time.5, method="kendall")
  constant.sims.skewness.5.tau[i] <- cor(constant.sims.skewness.5[,i], time.5, method="kendall")
  constant.sims.kurtosis.5.tau[i] <- cor(constant.sims.kurtosis.5[,i], time.5, method="kendall")
  constant.sims.id.5.tau[i] <- cor(constant.sims.id.5[,i], time.5, method="kendall")
  constant.sims.cor.5.tau[i] <- cor(constant.sims.cor.5[,i], time.5, method="kendall")
  constant.sims.cov.5.tau[i] <- cor(constant.sims.cov.5[,i], time.5, method="kendall")
  constant.sims.first.diff.variance.5.tau[i] <- cor(constant.sims.first.diff.variance.5[,i], time.5, method="kendall")
  constant.sims.decay.time.5.tau[i] <- cor(constant.sims.decay.time.5[,i], time.5, method="kendall")
  
}


# Compare Kendall Tau distributions #
par(mfrow=c(1,1))
hist(warming.sims.mean.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.mean.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cv.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cv.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.variance.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.variance.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.skewness.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.skewness.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.kurtosis.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.kurtosis.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.id.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.id.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cor.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cor.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cov.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cov.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.decay.time.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.decay.time.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.first.diff.variance.5.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.first.diff.variance.5.tau, col=alpha("blue", 0.4), breaks = 20,add=T)




#### From Brett et al.: The AUC can be efficiently calculated after ranking the combined set of test and nmull correlation coefficients by value,

# AUC = [rtest - ntest(nest+1)/2] / (ntest*nnull)
# where rtest is the sum of the ranks of test coefficients and 
# ntest and nnull are the number of realizations of the test and null models respectively

# Here is a function to calculate AUC based on 2 inputs:
# test.taus is a vector of Kendall tau values for the metric under warming/emerging conditions
# null.taus is a vector of Kendall tau values for the metric under constant/non-emerging conditions
# this will remove simulations that gave NAs for the Kendall tau correlation due to NAs in the metric

calc_auc <- function(test.taus, null.taus){
  
  combined.taus <- c(test.taus,null.taus)
  r <- rank(combined.taus, na.last=NA)
  
  test.taus.length <- sum(!is.na(test.taus))
  null.taus.length <- sum(!is.na(null.taus))
  
  r1 <- sum(r[1:test.taus.length])
  
  n1 <- test.taus.length
  n2 <- null.taus.length
  
  (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
}


AUC_variance_5 <- calc_auc(warming.sims.variance.5.tau,constant.sims.variance.5.tau) ; AUC_variance_5
AUC_mean_5 <- calc_auc(warming.sims.mean.5.tau,constant.sims.mean.5.tau); AUC_mean_5
AUC_kurtosis_5 <- calc_auc(warming.sims.kurtosis.5.tau,constant.sims.kurtosis.5.tau); AUC_kurtosis_5
AUC_skewness_5 <- calc_auc(warming.sims.skewness.5.tau,constant.sims.skewness.5.tau); AUC_skewness_5
AUC_id_5 <- calc_auc(warming.sims.id.5.tau,constant.sims.id.5.tau); AUC_id_5
AUC_first.diff.variance_5 <- calc_auc(warming.sims.first.diff.variance.5.tau,constant.sims.first.diff.variance.5.tau); AUC_first.diff.variance_5
AUC_cor_5 <- calc_auc(warming.sims.cor.5.tau,constant.sims.cor.5.tau); AUC_cor_5
AUC_cov_5 <- calc_auc(warming.sims.cov.5.tau,constant.sims.cov.5.tau); AUC_cov_5
AUC_cv_5 <- calc_auc(warming.sims.cv.5.tau,constant.sims.cv.5.tau); AUC_cv_5
AUC_decay.time_5 <- calc_auc(warming.sims.decay.time.5.tau,constant.sims.decay.time.5.tau); AUC_decay.time_5






# Window size of 15 #

warming.sims.variance.15 <- warming.sims.variance.15[1:45,]
warming.sims.mean.15 <- warming.sims.mean.15[1:45,]
warming.sims.cv.15 <- warming.sims.cv.15[1:45,]
warming.sims.skewness.15 <-  warming.sims.skewness.15[1:45,]
warming.sims.kurtosis.15 <-  warming.sims.kurtosis.15[1:45,]
warming.sims.id.15 <-  warming.sims.id.15[1:45,]
warming.sims.cor.15 <-  warming.sims.cor.15[1:45,]
warming.sims.cov.15 <-  warming.sims.cov.15[1:45,]
warming.sims.first.diff.variance.15 <-  warming.sims.first.diff.variance.15[1:45,]
warming.sims.decay.time.15 <-  warming.sims.decay.time.15[1:45,]


constant.sims.variance.15 <- constant.sims.variance.15[1:45,]
constant.sims.mean.15 <- constant.sims.mean.15[1:45,]
constant.sims.cv.15 <- constant.sims.cv.15[1:45,]
constant.sims.skewness.15 <-  constant.sims.skewness.15[1:45,]
constant.sims.kurtosis.15 <-  constant.sims.kurtosis.15[1:45,]
constant.sims.id.15 <-  constant.sims.id.15[1:45,]
constant.sims.cor.15 <-  constant.sims.cor.15[1:45,]
constant.sims.cov.15 <-  constant.sims.cov.15[1:45,]
constant.sims.first.diff.variance.15 <-  constant.sims.first.diff.variance.15[1:45,]
constant.sims.decay.time.15 <-  constant.sims.decay.time.15[1:45,]



### Time series for window of 5 ###
time.15 <- c(1:45)

warming.sims.variance.15.tau <- c()
warming.sims.mean.15.tau <- c()
warming.sims.cv.15.tau <- c()
warming.sims.skewness.15.tau <-  c()
warming.sims.kurtosis.15.tau <- c()
warming.sims.id.15.tau <-  c()
warming.sims.cor.15.tau <-  c()
warming.sims.cov.15.tau <-  c()
warming.sims.first.diff.variance.15.tau <-  c()
warming.sims.decay.time.15.tau <-  c()


constant.sims.variance.15.tau <- c()
constant.sims.mean.15.tau <- c()
constant.sims.cv.15.tau <- c()
constant.sims.skewness.15.tau <-  c()
constant.sims.kurtosis.15.tau <- c()
constant.sims.id.15.tau <- c()
constant.sims.cor.15.tau <- c()
constant.sims.cov.15.tau <- c()
constant.sims.first.diff.variance.15.tau <-  c()
constant.sims.decay.time.15.tau <- c()



# Calculate Kendall's tau for each metric vs the time series #
for(i in 1:num.of.sims){
  warming.sims.variance.15.tau[i] <- cor(warming.sims.variance.15[,i], time.15, method="kendall")
  warming.sims.mean.15.tau[i] <- cor(warming.sims.mean.15[,i], time.15, method="kendall")
  warming.sims.cv.15.tau[i] <- cor(warming.sims.cv.15[,i], time.15, method="kendall")
  warming.sims.skewness.15.tau[i] <- cor(warming.sims.skewness.15[,i], time.15, method="kendall")
  warming.sims.kurtosis.15.tau[i] <- cor(warming.sims.kurtosis.15[,i], time.15, method="kendall")
  warming.sims.id.15.tau[i] <- cor(warming.sims.id.15[,i], time.15, method="kendall")
  warming.sims.cor.15.tau[i] <- cor(warming.sims.cor.15[,i], time.15, method="kendall")
  warming.sims.cov.15.tau[i] <- cor(warming.sims.cov.15[,i], time.15, method="kendall")
  warming.sims.first.diff.variance.15.tau[i] <- cor(warming.sims.first.diff.variance.15[,i], time.15, method="kendall")
  warming.sims.decay.time.15.tau[i] <- cor(warming.sims.decay.time.15[,i], time.15, method="kendall")
  
  constant.sims.variance.15.tau[i] <- cor(constant.sims.variance.15[,i], time.15, method="kendall")
  constant.sims.mean.15.tau[i] <- cor(constant.sims.mean.15[,i], time.15, method="kendall")
  constant.sims.cv.15.tau[i] <- cor(constant.sims.cv.15[,i], time.15, method="kendall")
  constant.sims.skewness.15.tau[i] <- cor(constant.sims.skewness.15[,i], time.15, method="kendall")
  constant.sims.kurtosis.15.tau[i] <- cor(constant.sims.kurtosis.15[,i], time.15, method="kendall")
  constant.sims.id.15.tau[i] <- cor(constant.sims.id.15[,i], time.15, method="kendall")
  constant.sims.cor.15.tau[i] <- cor(constant.sims.cor.15[,i], time.15, method="kendall")
  constant.sims.cov.15.tau[i] <- cor(constant.sims.cov.15[,i], time.15, method="kendall")
  constant.sims.first.diff.variance.15.tau[i] <- cor(constant.sims.first.diff.variance.15[,i], time.15, method="kendall")
  constant.sims.decay.time.15.tau[i] <- cor(constant.sims.decay.time.15[,i], time.15, method="kendall")
  
}


# Compare Kendall Tau distributions #
par(mfrow=c(1,1))
hist(warming.sims.mean.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.mean.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cv.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cv.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.variance.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.variance.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.skewness.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.skewness.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.kurtosis.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.kurtosis.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.id.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.id.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cor.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cor.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cov.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cov.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.decay.time.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.decay.time.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.first.diff.variance.15.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.first.diff.variance.15.tau, col=alpha("blue", 0.4), breaks = 20,add=T)





AUC_variance_15 <- calc_auc(warming.sims.variance.15.tau,constant.sims.variance.15.tau) ; AUC_variance_15
AUC_mean_15 <- calc_auc(warming.sims.mean.15.tau,constant.sims.mean.15.tau); AUC_mean_15
AUC_kurtosis_15 <- calc_auc(warming.sims.kurtosis.15.tau,constant.sims.kurtosis.15.tau); AUC_kurtosis_15
AUC_skewness_15 <- calc_auc(warming.sims.skewness.15.tau,constant.sims.skewness.15.tau); AUC_skewness_15
AUC_id_15 <- calc_auc(warming.sims.id.15.tau,constant.sims.id.15.tau); AUC_id_15
AUC_first.diff.variance_15 <- calc_auc(warming.sims.first.diff.variance.15.tau,constant.sims.first.diff.variance.15.tau); AUC_first.diff.variance_15
AUC_cor_15 <- calc_auc(warming.sims.cor.15.tau,constant.sims.cor.15.tau); AUC_cor_15
AUC_cov_15 <- calc_auc(warming.sims.cov.15.tau,constant.sims.cov.15.tau); AUC_cov_15
AUC_cv_15 <- calc_auc(warming.sims.cv.15.tau,constant.sims.cv.15.tau); AUC_cv_15
AUC_decay.time_15 <- calc_auc(warming.sims.decay.time.15.tau,constant.sims.decay.time.15.tau); AUC_decay.time_15









# Window size of 30 #

warming.sims.variance.30 <- warming.sims.variance.30[1:30,]
warming.sims.mean.30 <- warming.sims.mean.30[1:30,]
warming.sims.cv.30 <- warming.sims.cv.30[1:30,]
warming.sims.skewness.30 <-  warming.sims.skewness.30[1:30,]
warming.sims.kurtosis.30 <-  warming.sims.kurtosis.30[1:30,]
warming.sims.id.30 <-  warming.sims.id.30[1:30,]
warming.sims.cor.30 <-  warming.sims.cor.30[1:30,]
warming.sims.cov.30 <-  warming.sims.cov.30[1:30,]
warming.sims.first.diff.variance.30 <-  warming.sims.first.diff.variance.30[1:30,]
warming.sims.decay.time.30 <-  warming.sims.decay.time.30[1:30,]


constant.sims.variance.30 <- constant.sims.variance.30[1:30,]
constant.sims.mean.30 <- constant.sims.mean.30[1:30,]
constant.sims.cv.30 <- constant.sims.cv.30[1:30,]
constant.sims.skewness.30 <-  constant.sims.skewness.30[1:30,]
constant.sims.kurtosis.30 <-  constant.sims.kurtosis.30[1:30,]
constant.sims.id.30 <-  constant.sims.id.30[1:30,]
constant.sims.cor.30 <-  constant.sims.cor.30[1:30,]
constant.sims.cov.30 <-  constant.sims.cov.30[1:30,]
constant.sims.first.diff.variance.30 <-  constant.sims.first.diff.variance.30[1:30,]
constant.sims.decay.time.30 <-  constant.sims.decay.time.30[1:30,]




### Time series for window of 30 ###
time.30 <- c(1:30)

warming.sims.variance.30.tau <- c()
warming.sims.mean.30.tau <- c()
warming.sims.cv.30.tau <- c()
warming.sims.skewness.30.tau <-  c()
warming.sims.kurtosis.30.tau <- c()
warming.sims.id.30.tau <-  c()
warming.sims.cor.30.tau <-  c()
warming.sims.cov.30.tau <-  c()
warming.sims.first.diff.variance.30.tau <-  c()
warming.sims.decay.time.30.tau <-  c()


constant.sims.variance.30.tau <- c()
constant.sims.mean.30.tau <- c()
constant.sims.cv.30.tau <- c()
constant.sims.skewness.30.tau <-  c()
constant.sims.kurtosis.30.tau <- c()
constant.sims.id.30.tau <- c()
constant.sims.cor.30.tau <- c()
constant.sims.cov.30.tau <- c()
constant.sims.first.diff.variance.30.tau <-  c()
constant.sims.decay.time.30.tau <- c()



# Calculate Kendall's tau for each metric vs the time series #
for(i in 1:num.of.sims){
  warming.sims.variance.30.tau[i] <- cor(warming.sims.variance.30[,i], time.30, method="kendall")
  warming.sims.mean.30.tau[i] <- cor(warming.sims.mean.30[,i], time.30, method="kendall")
  warming.sims.cv.30.tau[i] <- cor(warming.sims.cv.30[,i], time.30, method="kendall")
  warming.sims.skewness.30.tau[i] <- cor(warming.sims.skewness.30[,i], time.30, method="kendall")
  warming.sims.kurtosis.30.tau[i] <- cor(warming.sims.kurtosis.30[,i], time.30, method="kendall")
  warming.sims.id.30.tau[i] <- cor(warming.sims.id.30[,i], time.30, method="kendall")
  warming.sims.cor.30.tau[i] <- cor(warming.sims.cor.30[,i], time.30, method="kendall")
  warming.sims.cov.30.tau[i] <- cor(warming.sims.cov.30[,i], time.30, method="kendall")
  warming.sims.first.diff.variance.30.tau[i] <- cor(warming.sims.first.diff.variance.30[,i], time.30, method="kendall")
  warming.sims.decay.time.30.tau[i] <- cor(warming.sims.decay.time.30[,i], time.30, method="kendall")
  
  constant.sims.variance.30.tau[i] <- cor(constant.sims.variance.30[,i], time.30, method="kendall")
  constant.sims.mean.30.tau[i] <- cor(constant.sims.mean.30[,i], time.30, method="kendall")
  constant.sims.cv.30.tau[i] <- cor(constant.sims.cv.30[,i], time.30, method="kendall")
  constant.sims.skewness.30.tau[i] <- cor(constant.sims.skewness.30[,i], time.30, method="kendall")
  constant.sims.kurtosis.30.tau[i] <- cor(constant.sims.kurtosis.30[,i], time.30, method="kendall")
  constant.sims.id.30.tau[i] <- cor(constant.sims.id.30[,i], time.30, method="kendall")
  constant.sims.cor.30.tau[i] <- cor(constant.sims.cor.30[,i], time.30, method="kendall")
  constant.sims.cov.30.tau[i] <- cor(constant.sims.cov.30[,i], time.30, method="kendall")
  constant.sims.first.diff.variance.30.tau[i] <- cor(constant.sims.first.diff.variance.30[,i], time.30, method="kendall")
  constant.sims.decay.time.30.tau[i] <- cor(constant.sims.decay.time.30[,i], time.30, method="kendall")
  
}


# Compare Kendall Tau distributions #
par(mfrow=c(1,1))
hist(warming.sims.mean.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.mean.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cv.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cv.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.variance.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.variance.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.skewness.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.skewness.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.kurtosis.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.kurtosis.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.id.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.id.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cor.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cor.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.cov.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.cov.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.decay.time.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.decay.time.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)

hist(warming.sims.first.diff.variance.30.tau, col=alpha("red",0.4), breaks=20)
hist(constant.sims.first.diff.variance.30.tau, col=alpha("blue", 0.4), breaks = 20,add=T)





calc_auc <- function(test.taus, null.taus){
  
  combined.taus <- c(test.taus,null.taus)
  r <- rank(combined.taus, na.last=NA)
  
  test.taus.length <- sum(!is.na(test.taus))
  null.taus.length <- sum(!is.na(null.taus))
  
  r1 <- sum(r[1:test.taus.length])
  
  n1 <- test.taus.length
  n2 <- null.taus.length
  
  (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
}


AUC_variance_30 <- calc_auc(warming.sims.variance.30.tau,constant.sims.variance.30.tau) ; AUC_variance_30
AUC_mean_30 <- calc_auc(warming.sims.mean.30.tau,constant.sims.mean.30.tau); AUC_mean_30
AUC_kurtosis_30 <- calc_auc(warming.sims.kurtosis.30.tau,constant.sims.kurtosis.30.tau); AUC_kurtosis_30
AUC_skewness_30 <- calc_auc(warming.sims.skewness.30.tau,constant.sims.skewness.30.tau); AUC_skewness_30
AUC_id_30 <- calc_auc(warming.sims.id.30.tau,constant.sims.id.30.tau); AUC_id_30
AUC_first.diff.variance_30 <- calc_auc(warming.sims.first.diff.variance.30.tau,constant.sims.first.diff.variance.30.tau); AUC_first.diff.variance_30
AUC_cor_30 <- calc_auc(warming.sims.cor.30.tau,constant.sims.cor.30.tau); AUC_cor_30
AUC_cov_30 <- calc_auc(warming.sims.cov.30.tau,constant.sims.cov.30.tau); AUC_cov_30
AUC_cv_30 <- calc_auc(warming.sims.cv.30.tau,constant.sims.cv.30.tau); AUC_cv_30
AUC_decay.time_30 <- calc_auc(warming.sims.decay.time.30.tau,constant.sims.decay.time.30.tau); AUC_decay.time_30




### Summarize all of the AUC results in a data frame ###

metrics <- c("Variance","Mean","Kurtosis","Skewness","Ind. Disp.","First-Diff. Var.","AC","COV","CV","Decay Time")

AUC_5 <- c(AUC_variance_5,AUC_mean_5,AUC_kurtosis_5,AUC_skewness_5,AUC_id_5,AUC_first.diff.variance_5,AUC_cor_5,AUC_cov_5,AUC_cv_5,AUC_decay.time_5)
AUC_15 <- c(AUC_variance_15,AUC_mean_15,AUC_kurtosis_15,AUC_skewness_15,AUC_id_15,AUC_first.diff.variance_15,AUC_cor_15,AUC_cov_15,AUC_cv_15,AUC_decay.time_15)
AUC_30 <- c(AUC_variance_30,AUC_mean_30,AUC_kurtosis_30,AUC_skewness_30,AUC_id_30,AUC_first.diff.variance_30,AUC_cor_30,AUC_cov_30,AUC_cv_30,AUC_decay.time_30)


AUC_df <- data.frame(metrics,AUC_5,AUC_15,AUC_30); AUC_df








require(scales)
require(evobiR)
require(gdata)
require(doBy)
require(ecp)
require(strucchange)
require(TTR)
require(moments)


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


quartz()
par(mfrow=c(4,2))
par(oma=c(4,4,0,0))
par(mar=c(0,0,0,0))
plot(sub1$status.sum, type="l", col="blue", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
plot(sub2$status.sum,type="l", col="blue", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
plot(sub3$status.sum,type="l", col="blue", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
plot(sub4$status.sum,type="l", col="blue", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
plot(sub5$status.sum, type="l", col="red", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
plot(sub6$status.sum,type="l", col="red", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
plot(sub7$status.sum,type="l", col="red", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
plot(sub8$status.sum,type="l", col="red", xaxt="n", yaxt="n",ylim=c(0,6))
abline(v=20,lty=2)
mtext(text="Number of Infecteds",side=2, line=2, outer=TRUE)
mtext(text="Time (days)",side=1, line=2, outer=TRUE)



# Window length 5 (= 15 days) #
# Can't really do shorter window for our actual data, and can't really do longer or there is almost no baseline
# at all before R0 crosses 1 #


pop1.sd <- SlidingWindow(FUN="sd",data=sub1$status.sum, window=5,step=1)
pop2.sd <- SlidingWindow(FUN="sd",data=sub2$status.sum, window=5,step=1)
pop3.sd <- SlidingWindow(FUN="sd",data=sub3$status.sum, window=5,step=1)
pop4.sd <- SlidingWindow(FUN="sd",data=sub4$status.sum, window=5,step=1)
pop5.sd <- SlidingWindow(FUN="sd",data=sub5$status.sum, window=5,step=1)
pop6.sd <- SlidingWindow(FUN="sd",data=sub6$status.sum, window=5,step=1)
pop7.sd <- SlidingWindow(FUN="sd",data=sub7$status.sum, window=5,step=1)
pop8.sd <- SlidingWindow(FUN="sd",data=sub8$status.sum, window=5,step=1)

# Variance #
pop1.var <- pop1.sd^2
pop2.var <- pop2.sd^2
pop3.var <- pop3.sd^2
pop4.var <- pop4.sd^2
pop5.var <- pop5.sd^2
pop6.var <- pop6.sd^2
pop7.var <- pop7.sd^2
pop8.var <- pop8.sd^2

# Mean #
pop1.mean <- SlidingWindow(FUN="mean",data=sub1$status.sum, window=5,step=1)
pop2.mean <- SlidingWindow(FUN="mean",data=sub2$status.sum, window=5,step=1)
pop3.mean <- SlidingWindow(FUN="mean",data=sub3$status.sum, window=5,step=1)
pop4.mean <- SlidingWindow(FUN="mean",data=sub4$status.sum, window=5,step=1)
pop5.mean <- SlidingWindow(FUN="mean",data=sub5$status.sum, window=5,step=1)
pop6.mean <- SlidingWindow(FUN="mean",data=sub6$status.sum, window=5,step=1)
pop7.mean <- SlidingWindow(FUN="mean",data=sub7$status.sum, window=5,step=1)
pop8.mean <- SlidingWindow(FUN="mean",data=sub8$status.sum, window=5,step=1)

# Coefficient of variation, add 0.0000001 to the mean so that we don't divide by 0 #
pop1.cv <- pop1.sd/(pop1.mean+0.0000001)
pop2.cv <- pop2.sd/(pop2.mean+0.0000001)
pop3.cv <- pop3.sd/(pop3.mean+0.0000001)
pop4.cv <- pop4.sd/(pop4.mean+0.0000001)
pop5.cv <- pop5.sd/(pop5.mean+0.0000001)
pop6.cv <- pop6.sd/(pop6.mean+0.0000001)
pop7.cv <- pop7.sd/(pop7.mean+0.0000001)
pop8.cv <- pop8.sd/(pop8.mean+0.0000001)


# INDEX OF DISPERSION #
# index of dispersion = variance / mean #
pop1.id <- pop1.var/(pop1.mean+0.0000001)
pop2.id <- pop2.var/(pop2.mean+0.0000001)
pop3.id <- pop3.var/(pop3.mean+0.0000001)
pop4.id <- pop4.var/(pop4.mean+0.0000001)
pop5.id <- pop5.var/(pop5.mean+0.0000001)
pop6.id <- pop6.var/(pop6.mean+0.0000001)
pop7.id <- pop7.var/(pop7.mean+0.0000001)
pop8.id <- pop8.var/(pop8.mean+0.0000001)

### First-differenced variance ###
pop1.first.diff.variance <- c()
pop2.first.diff.variance <- c()
pop3.first.diff.variance <- c()
pop4.first.diff.variance <- c()
pop5.first.diff.variance <- c()
pop6.first.diff.variance <- c()
pop7.first.diff.variance <- c()
pop8.first.diff.variance <- c()


for(i in 2:35){
  pop1.first.diff.variance[i-1] <- pop1.var[i]-pop1.var[i-1]
  pop2.first.diff.variance[i-1] <- pop2.var[i]-pop2.var[i-1]
  pop3.first.diff.variance[i-1] <- pop3.var[i]-pop3.var[i-1]
  pop4.first.diff.variance[i-1] <- pop4.var[i]-pop4.var[i-1]
  pop5.first.diff.variance[i-1] <- pop5.var[i]-pop5.var[i-1]
  pop6.first.diff.variance[i-1] <- pop6.var[i]-pop6.var[i-1]
  pop7.first.diff.variance[i-1] <- pop7.var[i]-pop7.var[i-1]
  pop8.first.diff.variance[i-1] <- pop8.var[i]-pop8.var[i-1]
  
}

# Kurtosis #
pop1.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub1$status.sum, window=5,step=1)
pop2.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub2$status.sum, window=5,step=1)
pop3.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub3$status.sum, window=5,step=1)
pop4.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub4$status.sum, window=5,step=1)
pop5.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub5$status.sum, window=5,step=1)
pop6.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub6$status.sum, window=5,step=1)
pop7.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub7$status.sum, window=5,step=1)
pop8.kurtosis <- SlidingWindow(FUN="kurtosis",data=sub8$status.sum, window=5,step=1)

# Skewness #
pop1.skewness <- SlidingWindow(FUN="skewness",data=sub1$status.sum, window=5,step=1)
pop2.skewness <- SlidingWindow(FUN="skewness",data=sub2$status.sum, window=5,step=1)
pop3.skewness <- SlidingWindow(FUN="skewness",data=sub3$status.sum, window=5,step=1)
pop4.skewness <- SlidingWindow(FUN="skewness",data=sub4$status.sum, window=5,step=1)
pop5.skewness <- SlidingWindow(FUN="skewness",data=sub5$status.sum, window=5,step=1)
pop6.skewness <- SlidingWindow(FUN="skewness",data=sub6$status.sum, window=5,step=1)
pop7.skewness <- SlidingWindow(FUN="skewness",data=sub7$status.sum, window=5,step=1)
pop8.skewness <- SlidingWindow(FUN="skewness",data=sub8$status.sum, window=5,step=1)


# Now look for ACF-1 temporal autocorrelation #
# And Autocovariance #

pop1.ac <- c()
pop1.cov <- c()


for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub1$status.sum[c(window.start:window.end)]
  
  pop1.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop1.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
  
}


pop2.ac <- c()
pop2.cov <- c()

for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub2$status.sum[c(window.start:window.end)]
  
  pop2.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop2.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
}

pop3.ac <- c()
pop3.cov <- c()

for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub3$status.sum[c(window.start:window.end)]
  
  pop3.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop3.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
}


pop4.ac <- c()
pop4.cov <- c()

for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub4$status.sum[c(window.start:window.end)]
  
  pop4.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop4.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
}

pop5.ac <- c()
pop5.cov <- c()

for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub5$status.sum[c(window.start:window.end)]
  
  pop5.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop5.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
}

pop6.ac <- c()
pop6.cov <- c()

for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub6$status.sum[c(window.start:window.end)]
  
  pop6.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop6.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
}


pop7.ac <- c()
pop7.cov <- c()

for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub7$status.sum[c(window.start:window.end)]
  
  pop7.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop7.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
}

pop8.ac <- c()
pop8.cov <- c()

for(i in 1:35){
  window.start <- i
  window.end <- i+4
  
  segment <- sub8$status.sum[c(window.start:window.end)]
  
  pop8.ac[i] <- acf(segment, plot=FALSE, type="correlation")[[1]][2]
  pop8.cov[i] <- acf(segment, plot=FALSE, type="covariance")[[1]][2]
}



### Decay Time ###
# -time step / ln(autocorrelation) #
# This will create some NaNs, so use supress warnings function to surpress them #

pop1.decay.time <- c()
pop2.decay.time <- c()
pop3.decay.time <- c()
pop4.decay.time <- c()
pop5.decay.time <- c()
pop6.decay.time <- c()
pop7.decay.time <- c()
pop8.decay.time <- c()


for(i in 1:35){
  suppressWarnings(pop1.decay.time[i] <- -i/log(pop1.ac[i]))
  suppressWarnings(pop2.decay.time[i] <- -i/log(pop2.ac[i]))
  suppressWarnings(pop3.decay.time[i] <- -i/log(pop3.ac[i]))
  suppressWarnings(pop4.decay.time[i] <- -i/log(pop4.ac[i]))
  suppressWarnings(pop5.decay.time[i] <- -i/log(pop5.ac[i]))
  suppressWarnings(pop6.decay.time[i] <- -i/log(pop6.ac[i]))
  suppressWarnings(pop7.decay.time[i] <- -i/log(pop7.ac[i]))
  suppressWarnings(pop8.decay.time[i] <- -i/log(pop8.ac[i]))
}





### Get Kendall Tau statistics for each ###


### Time series for window of 5 ###
time.15 <- seq(1,15,1)

warming.variance.tau<- c(cor(pop5.var[1:15], time.15, method="kendall"),cor(pop6.var[1:15], time.15, method="kendall"),
                         cor(pop7.var[1:15], time.15, method="kendall"),cor(pop8.var[1:15], time.15, method="kendall"))
constant.variance.tau<- c(cor(pop1.var[1:15], time.15, method="kendall"),cor(pop2.var[1:15], time.15, method="kendall"),
                          cor(pop3.var[1:15], time.15, method="kendall"),cor(pop4.var[1:15], time.15, method="kendall"))


warming.mean.tau <- c(cor(pop5.mean[1:15], time.15, method="kendall"),cor(pop6.mean[1:15], time.15, method="kendall"),
                      cor(pop7.mean[1:15], time.15, method="kendall"),cor(pop8.mean[1:15], time.15, method="kendall"))
constant.mean.tau<- c(cor(pop1.mean[1:15], time.15, method="kendall"),cor(pop2.mean[1:15], time.15, method="kendall"),
                      cor(pop3.mean[1:15], time.15, method="kendall"),cor(pop4.mean[1:15], time.15, method="kendall"))


warming.cv.tau <- c(cor(pop5.cv[1:15], time.15, method="kendall"),cor(pop6.cv[1:15], time.15, method="kendall"),
                    cor(pop7.cv[1:15], time.15, method="kendall"),cor(pop8.cv[1:15], time.15, method="kendall"))
constant.cv.tau<- c(cor(pop1.cv[1:15], time.15, method="kendall"),cor(pop2.cv[1:15], time.15, method="kendall"),
                    cor(pop3.cv[1:15], time.15, method="kendall"),cor(pop4.cv[1:15], time.15, method="kendall"))


warming.skewness.tau <- c(cor(pop5.skewness[1:15], time.15, method="kendall"),cor(pop6.skewness[1:15], time.15, method="kendall"),
                          cor(pop7.skewness[1:15], time.15, method="kendall"),cor(pop8.skewness[1:15], time.15, method="kendall"))
constant.skewness.tau<- c(cor(pop1.skewness[1:15], time.15, method="kendall"),cor(pop2.skewness[1:15], time.15, method="kendall"),
                          cor(pop3.skewness[1:15], time.15, method="kendall"),cor(pop4.skewness[1:15], time.15, method="kendall"))


warming.kurtosis.tau <- c(cor(pop5.kurtosis[1:15], time.15, method="kendall"),cor(pop6.kurtosis[1:15], time.15, method="kendall"),
                          cor(pop7.kurtosis[1:15], time.15, method="kendall"),cor(pop8.kurtosis[1:15], time.15, method="kendall"))
constant.kurtosis.tau<- c(cor(pop1.kurtosis[1:15], time.15, method="kendall"),cor(pop2.kurtosis[1:15], time.15, method="kendall"),
                          cor(pop3.kurtosis[1:15], time.15, method="kendall"),cor(pop4.kurtosis[1:15], time.15, method="kendall"))


warming.id.tau <- c(cor(pop5.id[1:15], time.15, method="kendall"),cor(pop6.id[1:15], time.15, method="kendall"),
                    cor(pop7.id[1:15], time.15, method="kendall"),cor(pop8.id[1:15], time.15, method="kendall"))
constant.id.tau<- c(cor(pop1.id[1:15], time.15, method="kendall"),cor(pop2.id[1:15], time.15, method="kendall"),
                    cor(pop3.id[1:15], time.15, method="kendall"),cor(pop4.id[1:15], time.15, method="kendall"))


warming.ac.tau <- c(cor(pop5.ac[1:15], time.15, method="kendall"),cor(pop6.ac[1:15], time.15, method="kendall"),
                    cor(pop7.ac[1:15], time.15, method="kendall"),cor(pop8.ac[1:15], time.15, method="kendall"))
constant.ac.tau<- c(cor(pop1.ac[1:15], time.15, method="kendall"),cor(pop2.ac[1:15], time.15, method="kendall"),
                    cor(pop3.ac[1:15], time.15, method="kendall"),cor(pop4.ac[1:15], time.15, method="kendall"))


warming.cov.tau <- c(cor(pop5.cov[1:15], time.15, method="kendall"),cor(pop6.cov[1:15], time.15, method="kendall"),
                     cor(pop7.cov[1:15], time.15, method="kendall"),cor(pop8.cov[1:15], time.15, method="kendall"))
constant.cov.tau<- c(cor(pop1.cov[1:15], time.15, method="kendall"),cor(pop2.cov[1:15], time.15, method="kendall"),
                     cor(pop3.cov[1:15], time.15, method="kendall"),cor(pop4.cov[1:15], time.15, method="kendall"))


warming.first.diff.variance.tau <- c(cor(pop5.first.diff.variance[1:15], time.15, method="kendall"),cor(pop6.first.diff.variance[1:15], time.15, method="kendall"),
                                     cor(pop7.first.diff.variance[1:15], time.15, method="kendall"),cor(pop8.first.diff.variance[1:15], time.15, method="kendall"))
constant.first.diff.variance.tau<- c(cor(pop1.first.diff.variance[1:15], time.15, method="kendall"),cor(pop2.first.diff.variance[1:15], time.15, method="kendall"),
                                     cor(pop3.first.diff.variance[1:15], time.15, method="kendall"),cor(pop4.first.diff.variance[1:15], time.15, method="kendall"))


warming.decay.time.tau <- c(cor(pop5.decay.time[1:15], time.15, method="kendall"),cor(pop6.decay.time[1:15], time.15, method="kendall"),
                            cor(pop7.decay.time[1:15], time.15, method="kendall"),cor(pop8.decay.time[1:15], time.15, method="kendall"))
constant.decay.time.tau<- c(cor(pop1.decay.time[1:15], time.15, method="kendall"),cor(pop2.decay.time[1:15], time.15, method="kendall"),
                            cor(pop3.decay.time[1:15], time.15, method="kendall"),cor(pop4.decay.time[1:15], time.15, method="kendall"))



#### From Brett et al.: The AUC can be efficiently calculated after ranking the combined set of test and nmull correlation coefficients by value,

# AUC = [rtest - ntest(nest+1)/2] / (ntest*nnull)
# where rtest is the sum of the ranks of test coefficients and 
# ntest and nnull are the number of realizations of the test and null models respectively

# Here is a function to calculate AUC based on 2 inputs:
# test.taus is a vector of Kendall tau values for the metric under warming/emerging conditions
# null.taus is a vector of Kendall tau values for the metric under constant/non-emerging conditions
# this will remove simulations that gave NAs for the Kendall tau correlation due to NAs in the metric


calc_auc <- function(test.taus, null.taus){
  
  combined.taus <- c(test.taus,null.taus)
  r <- rank(combined.taus, na.last=NA)
  
  test.taus.length <- sum(!is.na(test.taus))
  null.taus.length <- sum(!is.na(null.taus))
  
  r1 <- sum(r[1:test.taus.length])
  
  n1 <- test.taus.length
  n2 <- null.taus.length
  
  (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
}


AUC_variance_data <- calc_auc(warming.variance.tau,constant.variance.tau) ; AUC_variance_data
AUC_mean_data <- calc_auc(warming.mean.tau,constant.mean.tau) ; AUC_mean_data
AUC_kurtosis_data <- calc_auc(warming.kurtosis.tau,constant.kurtosis.tau) ; AUC_kurtosis_data
AUC_skewness_data <- calc_auc(warming.skewness.tau,constant.skewness.tau) ; AUC_skewness_data
AUC_id_data <- calc_auc(warming.id.tau,constant.id.tau) ; AUC_id_data
AUC_first.diff.variance_data <- calc_auc(warming.first.diff.variance.tau,constant.first.diff.variance.tau) ; AUC_first.diff.variance_data
AUC_ac_data <- calc_auc(warming.ac.tau,constant.ac.tau) ; AUC_ac_data
AUC_cov_data <- calc_auc(warming.cov.tau,constant.cov.tau) ; AUC_cov_data
AUC_cv_data <- calc_auc(warming.cv.tau,constant.cv.tau) ; AUC_cv_data
AUC_decay.time_data <- calc_auc(warming.decay.time.tau,constant.decay.time.tau) ; AUC_decay.time_data



metrics <- c("Variance","Mean","Kurtosis","Skewness","Ind. Disp.","First-Diff. Var.","AC","COV","CV","Decay Time")

AUC_data <- c(AUC_variance_data,AUC_mean_data,AUC_kurtosis_data,AUC_skewness_data,AUC_id_data,AUC_first.diff.variance_data,AUC_ac_data,AUC_cov_data,AUC_cv_data,AUC_decay.time_data)


AUC_df <- data.frame(metrics,AUC_data); AUC_df










####### GET MEAN VALUES FROM SIMS FOR PLOT ######
# Get mean CV for warming and constant #
warming.sims.mean.incidence <- c()
for(i in 1:120){
  warming.sims.mean.incidence[i] <- mean(as.numeric(warming.incidence.daily[i,]))
}

constant.sims.mean.incidence <- c()
for(i in 1:120){
  constant.sims.mean.incidence[i] <- mean(as.numeric(constant.incidence.daily[i,]))
}


warming.sims.mean.variance <- c()
for(i in 1:45){
  warming.sims.mean.variance[i] <- mean(as.numeric(warming.sims.variance.15[i,]))
}

constant.sims.mean.variance <- c()
for(i in 1:45){
  constant.sims.mean.variance[i] <- mean(as.numeric(constant.sims.variance.15[i,]))
}


warming.sims.mean.mean <- c()
for(i in 1:45){
  warming.sims.mean.mean[i] <- mean(as.numeric(warming.sims.mean.15[i,]))
}

constant.sims.mean.mean <- c()
for(i in 1:45){
  constant.sims.mean.mean[i] <- mean(as.numeric(constant.sims.mean.15[i,]))
}







### CHANGE AXIS ON LAST COUPLE ON RIGHT AND GET TEXT LABELS TO LINE UP ##



######## PLOTTING INCIDENCE, MEAN, AND VARIANCE FOR DIFFERENT POPULATIONS ###

col.vec <- rep(c("gray0","gray12","gray24","gray36"),2)
time.plot <- c(6:20) # for experiments
plot.time <- c(16:60) # for sims

par(mfrow=c(3,4))
par(mar=c(1,1,1,1))
par(oma=c(4,4,4,0))

plot(constant.incidence.daily[,1] ~ c(1:120), type="n",ylim=c(0,100),xlim=c(0,120),
     xlab="Time",ylab="Disease Incidence",cex.lab=1.5, lwd=3)
for(i in 1:1000){
  lines(constant.incidence.daily[,i] ~ c(1:120),col=alpha("black",0.01), lwd=1)
}
lines(constant.sims.mean.incidence, col="black", lwd=3)
mtext("Disease Incidence",side=2, cex=1, line = 3, outer=FALSE)
mtext("Constant",side=3, cex=1, line = 1, outer=FALSE)
text(10, 95, labels=("a)"))


plot(warming.incidence.daily[,1], col="black", type="n",ylim=c(0,100),xlim=c(0,120),
     xlab="Time",ylab=NA,cex.lab=1.5)
for(i in 1:1000){
  lines(warming.incidence.daily[,i], col=alpha("black",0.01), lwd=1)
}
lines(warming.sims.mean.incidence, col="black", lwd=3)
abline(v=60, lty=2, col="yellow", lwd=2)
mtext("Warming",side=3, cex=1, line = 1, outer=FALSE)
text(10, 95, labels=("b)"))



plot(sub1$status.sum ~ sub1$day, type="n", ylim=c(0,6), ylab=NA, xlab=NA)
for(i in 1:4){
  sub <- subset(data.sum, population == i)
  lines(sub$status.sum ~ sub$day, col=alpha(col.vec[i],0.6), lty=1,lwd=1)
}
mtext("Constant",side=3, cex=1, line = 1, outer=FALSE)
text(5, 5.7, labels=("c)"))



plot(sub5$status.sum ~ sub5$day, type="n", ylim=c(0,6), ylab=NA, xlab=NA)
for(i in 5:8){
  sub <- subset(data.sum, population == i)
  lines(sub$status.sum ~ sub$day, col=alpha(col.vec[i],0.6), lty=1,lwd=1)
}
abline(v=c(60), col="yellow2",lty=2, lwd=2)
mtext("Warming",side=3, cex=1, line = 1, outer=FALSE)
mtext("Experiment",side=3, cex=1, line = 2, outer=TRUE, at=0.75)
mtext("Simulations",side=3, cex=1, line = 2, outer=TRUE, at=0.25)
text(5, 5.7, labels=("d)"))



# MEAN #

plot(constant.sims.mean.15[,1] ~ plot.time, col="black", type="n",ylim=c(0,30),xlim=c(0,120),
     xlab="Time",ylab="Mean",cex.lab=1.5)
for(i in 1:1000){
  lines(constant.sims.mean.15[,i] ~ plot.time, col=alpha("black",0.05), lwd=1)
}
lines(constant.sims.mean.mean ~ plot.time, col="black", lwd=3)
mtext("Mean",side=2, cex=1, line = 3, outer=FALSE)
text(5, 28.5, labels=("e)"))


plot(warming.sims.mean.15[,1] ~ plot.time, col="black", type="n",ylim=c(0,30),xlim=c(0,120),
     xlab="Time",ylab=NA,cex.lab=1.5)
for(i in 1:1000){
  lines(warming.sims.mean.15[,i] ~ plot.time, col=alpha("black",0.05), lwd=1)
}
lines(warming.sims.mean.mean ~ plot.time, col="black", lwd=3)
abline(v=60, lty=2, col="yellow", lwd=2)
text(5, 28.5, labels=("f)"))


plot(pop1.mean[1:15] ~ time.plot, col=alpha(col.vec[1],0.6),xaxt="n", type="l", ylim=c(0,1), xlim=c(0,40), lwd=1)
lines(pop2.mean[1:15] ~ time.plot, col=alpha(col.vec[2],0.6),  lwd=1, ylim=c(0,4))
lines(pop3.mean[1:15] ~ time.plot, col=alpha(col.vec[3],0.6), lwd=1,ylim=c(0,4))
lines(pop4.mean[1:15] ~ time.plot, col=alpha(col.vec[4],0.6), lwd=1,ylim=c(0,4))
axis(side=1, at = c(0,20/3,40/3,60/3,80/3,100/3,120/3), labels=c(0,20,40,60,80,100,120))
text(5, 0.95, labels=("g)"))


plot(pop5.mean[1:15] ~ time.plot, col=alpha(col.vec[1],0.6),xaxt="n", lwd=1,type="l",ylim=c(0,1), xlim=c(0,40))
lines(pop6.mean[1:15] ~ time.plot,col=alpha(col.vec[2],0.6), lwd=1, ylim=c(0,4))
lines(pop7.mean[1:15] ~ time.plot,col=alpha(col.vec[3],0.6), lwd=1, ylim=c(0,4))
lines(pop8.mean[1:15] ~ time.plot, col=alpha(col.vec[4],0.6),lwd=1, ylim=c(0,4))
axis(side=1, at = c(0,20/3,40/3,60/3,80/3,100/3,120/3), labels=c(0,20,40,60,80,100,120))
abline(v=c(20), col="yellow2",lty=2, lwd=2)
text(5, 0.95, labels=("h)"))


# variance #

plot(constant.sims.variance.15[,1] ~ plot.time, col="black", type="n",ylim=c(0,20),xlim=c(0,120),
     xlab="Time",ylab="Variance",cex.lab=1.5)
for(i in 1:1000){
  lines(constant.sims.variance.15[,i] ~ plot.time, col=alpha("black",0.01), lwd=1)
}
lines(constant.sims.mean.variance ~ plot.time, col="black", lwd=3)
mtext("Variance",side=2, cex=1, line = 3, outer=FALSE)
text(5, 19, labels=("i)"))


plot(warming.sims.variance.15[,1] ~ plot.time, col="black", type="n",ylim=c(0,20),xlim=c(0,120),
     xlab="Time",ylab=NA,cex.lab=1.5)
for(i in 1:1000){
  lines(warming.sims.variance.15[,i] ~ plot.time, col=alpha("black",0.01), lwd=1)
}
lines(warming.sims.mean.variance ~ plot.time, col="black", lwd=3)
abline(v=60, lty=2, col="yellow", lwd=2)
text(5, 19, labels=("j)"))


plot(pop1.var[1:15] ~ time.plot, col=alpha(col.vec[1],0.6), xaxt="n",lwd=1,type="l",ylim=c(0,1.5), xlim=c(0,40))
lines(pop2.var[1:15] ~ time.plot, col=alpha(col.vec[2],0.6), lwd=1, ylim=c(0,4))
lines(pop3.var[1:15] ~ time.plot, col=alpha(col.vec[3],0.6), lwd=1, ylim=c(0,4))
lines(pop4.var[1:15] ~ time.plot, col=alpha(col.vec[4],0.6), lwd=1, ylim=c(0,4))
axis(side=1, at = c(0,20/3,40/3,60/3,80/3,100/3,120/3), labels=c(0,20,40,60,80,100,120))
text(5, 1.425, labels=("k)"))


plot(pop5.var[1:15] ~ time.plot, col=alpha(col.vec[1],0.6), xaxt="n", lwd=1,type="l",ylim=c(0,1.5), xlim=c(0,40))
lines(pop6.var[1:15] ~ time.plot, col=alpha(col.vec[2],0.6), lwd=1, ylim=c(0,4))
lines(pop7.var[1:15] ~ time.plot, col=alpha(col.vec[3],0.6), lwd=1, ylim=c(0,4))
lines(pop8.var[1:15] ~ time.plot, col=alpha(col.vec[4],0.6), lwd=1, ylim=c(0,4))
axis(side=1, at = c(0,20/3,40/3,60/3,80/3,100/3,120/3), labels=c(0,20,40,60,80,100,120))
abline(v=c(20), col="yellow2",lty=2, lwd=2)
text(5, 1.425, labels=("l)"))


mtext("Time (days)",side=1, cex=1, line = 2, outer=TRUE)












