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
SamplingData <- read.csv("Kirk_Data/Epidemics_Data_Dryad/Kirk_et_al_epidemics_sampling_data.csv")
head(SamplingData)


