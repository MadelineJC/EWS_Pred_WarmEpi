require(scales)
require(evobiR)
require(gdata)
require(doBy)
require(ecp)
require(strucchange)
require(TTR)
require(moments)


data <- read.csv("Kirk_Data/Epidemics_Data_Dryad/Kirk_et_al_epidemics_sampling_data.csv")

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








#### CHECK AUC VALUES IF WE EXCLUDE POPULATION 5, WHICH NEVER TOOK OFF TO AN EPIDEMIC ###

### Get Kendall Tau statistics for each ###


### Time series for window of 5 ###
time.15 <- seq(1,15,1)

warming.variance.tau_without5 <- c(cor(pop6.var[1:15], time.15, method="kendall"),
                                   cor(pop7.var[1:15], time.15, method="kendall"),cor(pop8.var[1:15], time.15, method="kendall"))


warming.mean.tau_without5 <- c(cor(pop6.mean[1:15], time.15, method="kendall"),
                               cor(pop7.mean[1:15], time.15, method="kendall"),cor(pop8.mean[1:15], time.15, method="kendall"))

warming.cv.tau_without5 <- c(cor(pop6.cv[1:15], time.15, method="kendall"),
                             cor(pop7.cv[1:15], time.15, method="kendall"),cor(pop8.cv[1:15], time.15, method="kendall"))

warming.skewness.tau_without5 <- c(cor(pop6.skewness[1:15], time.15, method="kendall"),
                                   cor(pop7.skewness[1:15], time.15, method="kendall"),cor(pop8.skewness[1:15], time.15, method="kendall"))


warming.kurtosis.tau_without5 <- c(cor(pop6.kurtosis[1:15], time.15, method="kendall"),
                                   cor(pop7.kurtosis[1:15], time.15, method="kendall"),cor(pop8.kurtosis[1:15], time.15, method="kendall"))


warming.id.tau_without5 <- c(cor(pop6.id[1:15], time.15, method="kendall"),
                             cor(pop7.id[1:15], time.15, method="kendall"),cor(pop8.id[1:15], time.15, method="kendall"))


warming.ac.tau_without5 <- c(cor(pop6.ac[1:15], time.15, method="kendall"),
                             cor(pop7.ac[1:15], time.15, method="kendall"),cor(pop8.ac[1:15], time.15, method="kendall"))


warming.cov.tau_without5 <- c(cor(pop6.cov[1:15], time.15, method="kendall"),
                              cor(pop7.cov[1:15], time.15, method="kendall"),cor(pop8.cov[1:15], time.15, method="kendall"))


warming.first.diff.variance.tau_without5 <- c(cor(pop6.first.diff.variance[1:15], time.15, method="kendall"),
                                              cor(pop7.first.diff.variance[1:15], time.15, method="kendall"),cor(pop8.first.diff.variance[1:15], time.15, method="kendall"))


warming.decay.time.tau_without5 <- c(cor(pop6.decay.time[1:15], time.15, method="kendall"),
                                     cor(pop7.decay.time[1:15], time.15, method="kendall"),cor(pop8.decay.time[1:15], time.15, method="kendall"))


AUC_variance_data_without5 <- calc_auc(warming.variance.tau_without5 ,constant.variance.tau) ; AUC_variance_data_without5 
AUC_mean_data_without5  <- calc_auc(warming.mean.tau_without5 ,constant.mean.tau) ; AUC_mean_data_without5 
AUC_kurtosis_data_without5  <- calc_auc(warming.kurtosis.tau_without5, constant.kurtosis.tau) ; AUC_kurtosis_data_without5 
AUC_skewness_data_without5  <- calc_auc(warming.skewness.tau_without5 ,constant.skewness.tau) ; AUC_skewness_data_without5 
AUC_id_data_without5  <- calc_auc(warming.id.tau_without5 ,constant.id.tau) ; AUC_id_data_without5 
AUC_first.diff.variance_data_without5  <- calc_auc(warming.first.diff.variance.tau_without5 ,constant.first.diff.variance.tau) ; AUC_first.diff.variance_data_without5 
AUC_ac_data_without5  <- calc_auc(warming.ac.tau_without5 ,constant.ac.tau) ; AUC_ac_data_without5 
AUC_cov_data_without5  <- calc_auc(warming.cov.tau_without5 ,constant.cov.tau) ; AUC_cov_data_without5 
AUC_cv_data_without5  <- calc_auc(warming.cv.tau_without5 ,constant.cv.tau) ; AUC_cv_data_without5 
AUC_decay.time_data_without5  <- calc_auc(warming.decay.time.tau_without5 ,constant.decay.time.tau) ; AUC_decay.time_data_without5 


AUC_data <- c(AUC_variance_data_without5 ,AUC_mean_data_without5 ,AUC_kurtosis_data_without5 ,AUC_skewness_data_without5,
              AUC_id_data_without5 ,AUC_first.diff.variance_data_without5 ,AUC_ac_data_without5 ,AUC_cov_data_without5,
              AUC_cv_data_without5 ,AUC_decay.time_data_without5 )


AUC_df_without5 <- data.frame(metrics,AUC_data); AUC_df_without5







######## PLOTTING INCIDENCE, MEAN, AND VARIANCE FOR DIFFERENT POPULATIONS ###

col.vec <- rep(c("gray0","gray12","gray24","gray36"),2)
time.plot <- c(6:20)

quartz()
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
par(oma=c(4,4,0,0))
plot(sub1$status.sum ~ sub1$day, type="n", ylim=c(0,6), ylab=NA, xlab=NA, cex.axis=0.8)
for(i in 1:4){
  sub <- subset(data.sum, population == i)
  lines(sub$status.sum ~ sub$day, col=alpha(col.vec[i],0.6), lty=1,lwd=3)
}
text(5, 5.5, labels=("(a)"))
mtext("Disease Incidence",side=2, cex=1, line = 3, outer=FALSE)


plot(sub5$status.sum ~ sub5$day, type="n", ylim=c(0,6), ylab=NA, xlab=NA, cex.axis=0.8)
for(i in 5:8){
  sub <- subset(data.sum, population == i)
  lines(sub$status.sum ~ sub$day, col=alpha(col.vec[i],0.6), lty=1,lwd=3)
}
abline(v=c(15,30,45,75,90,105), lty=3)
abline(v=c(60), col="yellow2",lty=3, lwd=4)
text(x = c(7.5,22.5,37.5,52.5,67.5,82.5,97.5,112.5),y = 5.7, labels = c("10.0C","10.5C","11.0C","11.5C","12.0C","12.5C","13.0C","13.5C"), cex=0.75)
text(5, 5, labels=("(b)"))

# MEAN #
plot(pop1.mean[1:15] ~ time.plot, col=alpha(col.vec[1],0.6), type="l", ylim=c(0,1.5), xlim=c(0,40), lwd=3)
lines(pop2.mean[1:15] ~ time.plot, col=alpha(col.vec[2],0.6),  lwd=3, ylim=c(0,4))
lines(pop3.mean[1:15] ~ time.plot, col=alpha(col.vec[3],0.6), lwd=3,ylim=c(0,4))
lines(pop4.mean[1:15] ~ time.plot, col=alpha(col.vec[4],0.6), lwd=3,ylim=c(0,4))
mtext("Mean",side=2, cex=1, line = 3, outer=FALSE)


plot(pop5.mean[1:15] ~ time.plot, col=alpha(col.vec[1],0.6), lwd=3,type="l",ylim=c(0,1.5), xlim=c(0,40))
lines(pop6.mean[1:15] ~ time.plot,col=alpha(col.vec[2],0.6), lwd=3, ylim=c(0,4))
lines(pop7.mean[1:15] ~ time.plot,col=alpha(col.vec[3],0.6), lwd=3, ylim=c(0,4))
lines(pop8.mean[1:15] ~ time.plot, col=alpha(col.vec[4],0.6),lwd=3, ylim=c(0,4))
abline(v=c(20), col="yellow2",lty=3, lwd=4)


# variance #
plot(pop1.var[1:15] ~ time.plot, col=alpha(col.vec[1],0.6), lwd=3,type="l",ylim=c(0,2), xlim=c(0,40))
lines(pop2.var[1:15] ~ time.plot, col=alpha(col.vec[2],0.6), lwd=3, ylim=c(0,4))
lines(pop3.var[1:15] ~ time.plot, col=alpha(col.vec[3],0.6), lwd=3, ylim=c(0,4))
lines(pop4.var[1:15] ~ time.plot, col=alpha(col.vec[4],0.6), lwd=3, ylim=c(0,4))
mtext("Variance",side=2, cex=1, line = 3, outer=FALSE)


plot(pop5.var[1:15] ~ time.plot, col=alpha(col.vec[1],0.6), lwd=3,type="l",ylim=c(0,2), xlim=c(0,40))
lines(pop6.var[1:15] ~ time.plot, col=alpha(col.vec[2],0.6), lwd=3, ylim=c(0,4))
lines(pop7.var[1:15] ~ time.plot, col=alpha(col.vec[3],0.6), lwd=3, ylim=c(0,4))
lines(pop8.var[1:15] ~ time.plot, col=alpha(col.vec[4],0.6), lwd=3, ylim=c(0,4))
abline(v=c(60), col="yellow2",lty=3, lwd=4)

mtext("Time (days)",side=1, cex=1, line = 2, outer=TRUE)






