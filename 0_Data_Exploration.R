#### Required packages ####
library(tidyverse)

#### Kirk_et_al_epidemics_sampling_data.csv ####

## About:
# •	population: population replicate number
# •	treatment: stable temperature or rising temperature
# •	temp: temperature in °C
# •	period: number (1-8) representing a fifteen day temperature period for warming treatment
# •	day: day of the experiment
# •	date: date
# •	rep: number for individual Daphnia sampled from each population each day
# •	status: infection status, 0 = uninfected, 1 = infected
# •	intensity: intensity of infection (0 for uninfected individuals)

## Import
sampling_data <- read.csv("Kirk_Data/Epidemics_Data_Dryad/Kirk_et_al_epidemics_sampling_data.csv")
head(sampling_data)
nrow(sampling_data)

## Explore
unique(sampling_data$population) # Eight replicate populations
unique(sampling_data$treatment) # Two temperature treatments; stable and rising
unique(sampling_data$temp) # Eight different temperatures; 10C to 13.5C in incrememnts of 0.5C
max(sampling_data$day) # 120 day experiment
length(unique(sampling_data$rep)) # 12 Daphnia samples from each population, per day
summary(sampling_data$intensity) # Min. = 0; Median = 0; Mean = 1.32; Max. = 183

constant <- subset(sampling_data, treatment == 'stable')
constant_pop1 <- subset(constant, population == 1)
constant_pop2 <- subset(constant, population == 2)
constant_pop3 <- subset(constant, population == 3)
constant_pop4 <- subset(constant, population == 4)

rising <- subset(sampling_data, treatment == 'rising')
rising_pop1 <- subset(rising, population == 1)
rising_pop2 <- subset(rising, population == 2)
rising_pop3 <- subset(rising, population == 3)
rising_pop4 <- subset(rising, population == 4)

## Some plots
colours1 <- c("navy", "royalblue", "lightblue", "lightblue1")
colours2 <- c("#9F2B68", "#DE3163", "#FF69B4", "#F8C8DC")

par(mfrow = c(1, 1))
plot(constant$day, constant$intensity, 
     pch = 16, col = adjustcolor(colours1[factor(constant$population)], alpha.f = 0.5),
     xlab = "Day", ylab = "Infection Burden", 
     ylim = c(0, 200))
points(rising$day, rising$intensity,
       pch = 17, col = adjustcolor(colours2[factor(rising$population)], alpha.f = 0.5))

par(mfrow = c(1, 2))
plot(constant$day, constant$intensity, 
     pch = 16, col = adjustcolor(colours1[factor(constant$population)], alpha.f = 0.7),
     main = "Constant Temperature", xlab = "Day", ylab = "Infection Burden",
     ylim = c(0, 200))
legend("topleft", legend = c("Population 1", "Population 2", "Population 3", "Population 4"),
       pch = 16, col = colours1, cex = 1, pt.cex = 2, bty = "n", x.intersp = 0.5)
plot(rising$day, rising$intensity, 
     pch = 17, col = adjustcolor(colours2[factor(rising$population)], alpha.f = 0.7),
     main = "Rising Temperature", xlab = "Day", ylab = "Infection Burden",
     ylim = c(0, 200))
legend("topleft", legend = c("Population 5", "Population 6", "Population 7", "Population 8"),
       pch = 16, col = colours2, cex = 1, pt.cex = 2, bty = "n", x.intersp = 0.5)

#### Plotting Simulated and Empirical Data ####
sim.N <- 170; samp.N <- 12; days <- seq(1, 120, 1); exp.days <- seq(3, 120, 3)
sim.c <- cbind(constant1[ , 2:101], constant2[ , 2:101], constant3[ , 2:101], constant4[ , 2:101], constant5[ , 2:101], constant6[ , 2:101], constant7[ , 2:101], constant8[ , 2:101], constant9[ , 2:101], constant10[ , 2:101])
sim.w <- cbind(warming1[ , 2:101], warming2[ , 2:101], warming3[ , 2:101], warming4[ , 2:101], warming5[ , 2:101], warming6[ , 2:101], warming7[ , 2:101], warming8[ , 2:101], warming9[ , 2:101], warming10[ , 2:101])
exp.c <- cbind(sub1$status, sub2$status, sub3$status, sub4$status)
exp.w <- cbind(sub5$status, sub6$status, sub7$status, sub8$status)
sim.c.prop <- sim.c/sim.N; sim.w.prop <- sim.w/sim.N
exp.c.prop <- exp.c/samp.N; exp.w.prop <- exp.w/samp.N

plot(1, type = "n", main = "Constant Treatment", xlab = "Day", ylab = "Prevalence", xlim = c(1, 120), ylim = c(0, 1))
for (i in 1:1000){
        lines(x = days, y = sim.c.prop[ , i], col = adjustcolor("royalblue", alpha = 0.1))
}
lines(x = days, y = rowMeans(sim.c.prop), col = "lightblue1", lwd = 3)
for (i in 1:4){
        lines(x = exp.days, y = exp.c.prop[ , i], col = "navy", lwd = 3)
}

plot(1, type = "n", main = "Warming Treatment", xlab = "Day", ylab = "Prevalence", xlim = c(1, 120), ylim = c(0, 1))
for (i in 1:1000){
        lines(x = days, y = sim.w.prop[ , i], col = adjustcolor("#DE3163", alpha = 0.1))
}
lines(x = days, y = rowMeans(sim.w.prop), col = "orchid1", lwd = 3)
for (i in 1:4){
        lines(x = exp.days, y = exp.w.prop[ , i], col = "#9F2B68", lwd = 3)
}
abline(v = 60, col = "orange1", lwd = 3)









