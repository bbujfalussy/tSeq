########################################################################
## read data from the matlab file

library(R.matlab)
source('./bandFilt.R', chdir = TRUE)

sampling.rate <- 1250 # Hz

i_rat <- 1
i_day <- 1

lfp_file <- paste('./Pfeiffer13/Rat', i_rat, '_Day', i_day, '_LFP_', lfp_side, '.mat', sep='')
tryCatch({lfp_data <- readMat(lfp_file)}, error = function(e) print('error: LFP file not found. Please contact the authors to obtain the data.'), warning=function(w) print('error: LFP file not found.  Please contact the authors to obtain the data.'))

summary(data)
data$LFP.Electrodes
data$LFP.Frequency
dim(data$LFP.Left)
freq <- as.vector(data$LFP.Frequency)
dt <- 1/freq # s

Tdur <- 60 #s
N <- Tdur * freq

lfp <- data$LFP.Left[1:(N+1),2]

lfp.theta <- bandFilt(lfp, dt, f.low=4, f.high=12, graphics=T, fast=F)
lfp.ripple <- bandFilt(lfp, dt, f.low=120, f.high=250, graphics=T, fast=F)
tt <- seq(0, Tdur, by=dt)


matplot(tt, cbind(lfp, lfp.theta, lfp.ripple), col=c(1,2,3), lty=1, t='l')

for (j in 1:100){
	tt <- j * 30000 + 1:30000
	plot(data$LFP.Left[tt,1], data$LFP.Left[tt,2], t='l')
	readline(j)
}


