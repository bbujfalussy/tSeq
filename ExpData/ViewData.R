########################################################################
## read data from the matlab file

library(R.matlab)

sampling.rate <- 1250 # Hz

tryCatch({data <- readMat("Pfeiffer13/spikes.mat")}, error = function(e) print('File containing the spikes not found. Please contact the authors to obtain it!'), warning=function(w) print('File containing the spikes not found. Please contact the authors to obtain it!'))
data <- data[[1]]
names(data) <- list('rat1', 'rat2', 'rat3', 'rat4')
summary(data)
for (i in 1:4){
	names(data[[i]]) <- list('day1', 'day2')
	for (j in 1:2) names(data[[i]][[j]]) <- list('spikes', 'position', 'I_Cells', 'E_Cells')
	
} 

### contains 4 lists, corresponding to 4 rats

### each rat has 2 experimental days

## for each day the data contains
## 1. matrix, SpikeData is a time-ordered list of spikes, with the first column being the timepoint (in seconds) of each spike and the second column being the cell ID.  1000 Hz
# 2. PositionData is a time-ordered list of positions and head directions, with the first column being the timepoint (in seconds, synchronized with the spike data) of each camera frame, the second column being the rat's X position (in cm), the third column being the rat's Y position (in cm), and the fourth column being the rat's head direction (in degrees).  30 Hz
# 3-4: ExcitatoryNeurons/InhibitoryNeurons is a list of the cell IDs for neurons that were classified as excitatory vs. inhibitory.

plot(data$rat1$day2$spikes[1:10000, 1], data$rat1$day2$spikes[1:10000, 2], pch=16, cex=0.3, xlab='time (s)', ylab='cells')
abline(h=data$rat1$day2$I_Cells, col=2)

plot(data$rat1$day1$position[1:1000, 1], data$rat1$day1$position[1:1000, 2], t='l')
plot(data$rat1$day1$position[, 1], data$rat1$day1$position[, 2], t='l')
plot(data$rat1$day1$position[, 1], data$rat1$day1$position[, 3], t='l')

plot(data$rat1$day1$position[, 2], data$rat1$day1$position[, 3], t='l')
lines(data$rat1$day2$position[, 2], data$rat1$day2$position[, 3], t='l', col=3)

plot(data$rat1$day2$position[, 2], data$rat1$day2$position[, 3], t='l', col=3)

