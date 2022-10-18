##############################################
## synthetic data for theta sequences in the hippocampus
##############################################
## for the generative model of the motion see kalman_motion.pdf
## STEPS I will need
## 1. generate data using smooth, continuous trajectories
## 2. inference runs for a given set of observations
## 3. sampling
## 4. recursive inference + sampling

## 1. Generate data
graphics <- F
set.seed(33)
library(Matrix)
library(viridis)
require(gplots)
source('~/Projects/PATTERN/ThetaSeq/gen_inf_kalman.R', chdir = TRUE)
source('~/Programs/Code/HPC_Data/SimRunND.R', chdir = TRUE)
source('~/Projects/PATTERN/ThetaSeq/FosterData/PlaceCellFunctions.R', chdir = TRUE)
########################################################
# 1. generating realistic trajectories - 3 Hz, smooth speed 
## to simulate motion, we will do the following steps - kalman_motion.pdf:
## 1.1. Simulate motion in 2D, using smooth velocity changes. This provides the intended position for each timestep - x.intended
## 1.2. Generate an initial true position x_1
## 1.3. Generate an observation y_1
## 1.4. Use filtering to calculate posterior mean of the position
## 1.5. Calculate the drive for the next step as the difference between the intended position and the posterior mean position 
## 1.6. The new position will be computed from old position, the external drive and the noise
## repeat 1.3 - 1.6

Tmax <- 60*60 # 5 min
dt.run <- 1/10 # s, resolution should be matched to theta frequency
T.axis <- seq(dt.run, Tmax, by= dt.run)

Nmax <- Tmax / dt.run
L <- 2 # m box size
alpha <- 0.25 # 0.25
# parameters of speed - the absolute value of velocity
pars.speed <- list(mu.v=0.16, # m/s (in 1D, animals move much faster, 20-40 cm/s is typical. Also rats are faster than mice.)
	sd.v=0.1, # m/s, sd of the speed
	tau.v=2, # 2 - s
	v.min=0.02, # 0.02 min of speed
	v.max=0.8, # max of speed
	m.dv=0.035, # 0.05 change in the direction of motion
	alpha=alpha, # smoothness of the drive
	dt=dt.run, # time step
	xmax=L) # arena size

q <- 0.015^2* dt.run/0.1 # evolution noise - 1.5 cm - should be smaller than the observation noise
r <- 0.15^2* 0.1 / dt.run # observation noise - 15 cm

ss.filt <- (sqrt(q^2 + 4*r*q) - q ) / 2
# c(sqrt(ss.filt), sqrt(min(filtered.post$var)))

ss.smooth <- ss.filt*(ss.filt+q) / (q + 2* ss.filt)
# c(sqrt(ss.smooth), sqrt(min(smooth.post$var)))

c(sqrt(ss.smooth), sqrt(ss.filt), sqrt(r))

########################################
## 1.1. Simulate motion in 2D, using smooth velocity changes. This is going to be the intended position x(t).
## vx$x : matrix, each column is a position coordinate
## vx$v : matrix, each column is a velocity
vx_raw <- sim.run(Tmax=Tmax, xmax=c(2, 0.1), dt=dt.run, pars.speed, ndim=2, seed=11) 

fxpos <- vx_raw$x[1,]
TT <- length(fxpos)
tpos <- seq(0, by=dt.run, length=TT)

th_right <- 1.98
th_left <- 0.02

istart.right <- which((fxpos[-TT] < th_left) & (fxpos[-1] > th_left))
iend.right <- which((fxpos[-TT] < th_right) & (fxpos[-1] > th_right))
istart.left <- which((fxpos[-TT] > th_right) & (fxpos[-1] < th_right))
iend.left <- which((fxpos[-TT] > th_left) & (fxpos[-1] < th_left))

right.runs <- matrix(c(0,0), 1)
ii_run <- 1
for (i.run in 1:length(istart.right)){
	istart <- istart.right[i.run]
	ii <- which(iend.right > istart)
	if (length(ii) > 0){
		iend <- iend.right[min(ii)]
		if ((iend - istart) <= 150){
			N <- nrow(right.runs) # check if it returned to start ...
			if (iend == right.runs[N,2]){
				right.runs[N,] <- c(istart, iend)
			} else {
				right.runs <- rbind(right.runs, c(istart, iend))
			}
		}
	}
}
right.runs <- right.runs[-1,]

left.runs <- matrix(c(0,0), 1)
for (i.run in 1:length(istart.left)){
	istart <- istart.left[i.run]
	ii <- which(iend.left > istart)
	if (length(ii) > 0){
		iend <- iend.left[min(ii)]
		if ((iend - istart) <= 150){
			N <- nrow(left.runs) # check if it returned to start ...
			if (iend == left.runs[N,2]){
				left.runs[N,] <- c(istart, iend)
			} else {
				left.runs <- rbind(left.runs, c(istart, iend))
			}
		}
	}
}
left.runs <- left.runs[-1,]

plot(tpos, fxpos, col=1, t='l', xlim=c(0, 200))
abline(v=tpos[right.runs[,1]], col=3)
abline(v=tpos[right.runs[,2]], col=2)

abline(v=tpos[left.runs[,1]], col=4, lty=2)
abline(v=tpos[left.runs[,2]], col=6, lty=2)


vx <- list(x=vx_raw$x[,right.runs[1,1]:right.runs[1,2]], v=vx_raw$v[,right.runs[1,1]:right.runs[1,2]])
imax <- right.runs[1,2]
add_more <- T
while(add_more){
	ii <- min(which(left.runs[,1] > imax))
	if (ii == Inf) break
	vx$x <- cbind(vx$x, vx_raw$x[,left.runs[ii,1]:left.runs[ii,2]])	
	vx$v <- cbind(vx$v, vx_raw$v[,left.runs[ii,1]:left.runs[ii,2]])	
	imax <- left.runs[ii,2]

	ii <- min(which(right.runs[,1] > imax))
	if (ii == Inf) break
	vx$x <- cbind(vx$x, vx_raw$x[,right.runs[ii,1]:right.runs[ii,2]])	
	vx$v <- cbind(vx$v, vx_raw$v[,right.runs[ii,1]:right.runs[ii,2]])	
	imax <- right.runs[ii,2]
}

plot(vx$x[1,1:500], vx$x[2,1:500], t='l')
Tmax <- ncol(vx$x) / 10
T.axis <- seq(dt.run, Tmax, by= dt.run)
Nmax <- round(Tmax / dt.run)
########################################
## 1.2-1.5. Simulate motion in 2D, using smooth velocity changes.
## now we simulate the Kalman step one by one
xy <- list(x=matrix(NA, 2, Nmax), y=matrix(NA, 2, Nmax))
filtered.post <- list(mu=matrix(NA, 2, Nmax), var=matrix(NA, 2, Nmax))
drive <- matrix(0, 2, Nmax)

# 1.2 -1.3
## sim.kalman.step(current position, evolution noise, observation noise, drive)
## we start with a zero drive...
xy.t <- sim.kalman.step(vx$x[,1], q, r, 0)
prior <- list(mu=vx$x[,1], var=c(r,r))
# 1.4
filtered.post.t <- filt.kalman(matrix(xy.t$y, 2), q, r, prior)

xy$x[,1] <- xy.t$x
xy$y[,1] <- xy.t$y
filtered.post$mu[,1] <- filtered.post.t$mu
filtered.post$var[,1] <- filtered.post.t$var
drive.t <- matrix(0, 2, 1)

## start from the true position, and perform the filtering and the position update sequentially
for (i.time in 2:Nmax){
	drive.t <- drive.t * (1-alpha) + alpha * (vx$x[,i.time] - filtered.post.t$mu)
	# 1.2 -1.3
	xy.t <- sim.kalman.step(xy.t$x, q, r, drive.t)
	# 1.4
	filtered.post.t <- filt.kalman(matrix(xy.t$y, 2), q, r, filtered.post.t, drive.t)

	## storing the variables	
	xy$x[,i.time] <- xy.t$x
	xy$y[,i.time] <- xy.t$y
	filtered.post$mu[,i.time] <- filtered.post.t$mu
	filtered.post$var[,i.time] <- filtered.post.t$var
	drive[,i.time] <- drive.t
	
} 


## now we run the local inference: at each time step we 
# 1. update the filtering 
# 2. calculate smoothing back to a couple of timesteps for t.past
# 3. make predictions into the future for t.future
# the local posteriors are represented by samples 
# IMPORTANT: when sim_irreg == T: 
# each cycle is different:
# - they start at different point in the past -1.5 < t.past < 0 
# - they move to a different distance into the future 0 < t.future < 1.5
# - with encoding 1 < t.total < 3 s 
# - each cycle lasts for 80 < t.theta < 160 ms - so the encoding speed is also variable - built in the local.infer.kalman function
#
# when sim_irreg == F:

sim.irreg <- F
N.samples <- 100
outfile <- paste('motiondata1D_Tmax', Tmax, '.RData', sep='')
if (outfile %in% list.files()){
	load(file=outfile)
} else {
	post <- local.infer.kalman(xy$y, q, r, drive, N.samples, pars.speed, X=xy$x, graphics=F, irreg=FALSE)
	# each theta cycle represents the inferred trajectory up to now
	# the last theta cycle starts when sensory data ends
	motiondata <- list(T.axis=T.axis, intended=vx, position=xy, post=post, pars.speed=pars.speed, q=q, r=r, drive=drive)
	save(motiondata, file=outfile)
}

plot(xy$x[1,1:2000], xy$x[2,1:2000], t='o', pch=16, col=viridis(2000))
plot(T.axis[1:2000], xy$x[1,1:2000], t='o', pch=16, col=viridis(2000))

##############################################
# # 1. no uncertainty is represented in the population activity (PA)
# # 		- PA represents the the mean
# # 		- the error between two subsequent trajectories is SMALLER than the error between the true and a single trajectory
# # 		- trajectories are necoded with CONSTANT precision
# #
# # 2. sampling
# #		- PA samples from the posterior
# # 		- the error between two subsequent trajectories is SIMILAR to the error between the true and one of the trajectories
# # 		- trajectories are necoded with CONSTANT precision
# #
# # # 3. PPC
# # 		- the PA encodes the full posterior at any given time
# # 		- the error between two subsequent trajectories is SMALLER than the error between the true and a single trajectory
# # 		- trajectories are necoded with DECREASING precision
# #		- variance is encoded in the firing rate
# #
# # # 4. DDC
# # 		- the PA encodes the full posterior at any given time
# # 		- the error between two subsequent trajectories is SMALLER than the error between the true and a single trajectory
# # 		- trajectories are necoded with DECREASING precision
# # 		- variance is encoded in the population activity

###############################################
## Generate synthetic neuronal data to encode theta sequences

## 1. generate place cell tuning
N.cells <- 200
cat("generating the rate matrix... \n")
# png(file='Figs/placecells16.png', 800, 800)
## Mixture of Gaussian place fields: more realistic, but inference is slower - can only be implemented using the histogram method
pcs.rate <- gen.rate(N.cells, 2, 1, 10, seed=19, active.cells=T, rate.max=c(5, 15),  rate0=c(0.1,0.25), r.place.field=c(0.1, 0.3), graphics=0, rand.nfields=F)
## Multivariate Gaussian place fields: inference can be done analytically
# pcs.rate <- gen.rate(N.cells, 2, 81, 1, seed=19, active.cells=T, rate.max=c(5, 30),  rate0=c(0.1,0.25), r.place.field=c(0.1, 0.344), graphics=F, rand.nfields=F)
# dev.off()

unsorted_rates <- array(NA, dim=c(N.cells, 44, 10))
xx <- rep(seq(-7.5, 207.5, by=5)/100, 10)
yy <- rep(seq(-12.5, 32.5, by=5)/100, each=44)
for (i in 1:N.cells) unsorted_rates[i,,] <- t(matrix(get.rate.x(xx=rbind(xx, yy), pcs.rate[[i]]), 10, byrow=T))
# hist(apply(rates, 1, mean), br=seq(0, 30, by=2))
# hist(apply(rates*0.1, c(2,3), sum))

i_sort <- sort(apply(unsorted_rates[,,5], 1, which.max), ind=T)$ix
rates <- unsorted_rates[i_sort,,]

matplot(t(rates[,,5]), t='l', lty=1, col=viridis(200))
image(1:44, 1:200, t(rates[,,5]), col=viridis(24))

spPth <- round(mean(apply(rates*0.1, c(2,3), sum))) # spikes per theta - should be 31
cat('spikes per theta: ', spPth, '±', sd(apply(rates*0.1, c(2,3), sum)))


attr(rates, 'xcenters') <- seq(-12.5, 32.5, by=5)
attr(rates, 'ycenters') <- seq(-7.5, 207.5, by=5)

## 1. generate place cell ratemaps
t.samples <- seq(1, 100, length=21)
Nfilt <- 24
filt <- rep(1, Nfilt)/Nfilt

i.code <- 5
codes <- c('MAP', 'sampling', 'PPC', 'DDC', 'PPC_var')
N.theta <- length(motiondata$post)

for (i.code in 1:4){
	code <- codes[i.code]
	# code <- 'DDC' # sampling, PPC, DDC, MAP
	reftau <- 30
	spt <- c(0, 0, 0) # spike times, cellids, spike phases
	xy.trajectories <- list() # the smooth trajectory, used for generating place cell activity
	var.trajectories <- list()
	
	t0 <- motiondata$post[[1]]$tstats[1]
	theta.starts <- t0

	for (i.theta in 1:N.theta) { # 191: time after which the smoothing achieved steady state?
	
		xy.tau <- motiondata$post[[i.theta]]$estimpar$mu # for MAP, DDCal and PPC, the posterior in "motion-coordinates"
		var.xy.tau <- motiondata$post[[i.theta]]$estimpar$var # for PPC		
		if (code=='sampling')	{
			xy.tau <- motiondata$post[[i.theta]]$samples[1,,] # for sampling	
		} 
		L.theta <- ncol(xy.tau)
		
		T.theta <- motiondata$post[[i.theta]]$tstats[4] * 1000 # this is the real duration of the theta sequence in ms
		t.samples <- seq(1, T.theta, length=L.theta)
		t.samples.ms <- seq(1, T.theta)
		
		xy.ms <- rbind(approx(t.samples, xy.tau[1,], t.samples.ms)$y, approx(t.samples, xy.tau[2,], t.samples.ms)$y) # interpolating
		xy.ms <- cbind(matrix(rep(xy.ms[,1], Nfilt/2), 2), xy.ms, matrix(rep(xy.ms[,T.theta], Nfilt/2), 2)) # padding for filtering
		xy.ms <- rbind(filter(xy.ms[1,], filt), filter(xy.ms[2,], filt))[,(Nfilt/2):(T.theta+Nfilt/2-1)] # filtering
		# plot(xy.ms[1,], xy.ms[2,])
		xy.trajectories[[i.theta]] <- xy.ms
		# plot(xy.ms[1,], xy.ms[2,], xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), t='l')
	
		var.xy.ms <- rbind(approx(t.samples, var.xy.tau[1,], t.samples.ms)$y, approx(t.samples, var.xy.tau[2,], t.samples.ms)$y) # interpolating
		var.xy.ms <- cbind(matrix(rep(var.xy.ms[,1], Nfilt/2), 2), var.xy.ms, matrix(rep(var.xy.ms[, T.theta], Nfilt/2), 2)) # padding for filtering
		var.xy.ms <- rbind(filter(var.xy.ms[1,], filt), filter(var.xy.ms[2,], filt))[,(Nfilt/2):(T.theta +Nfilt/2-1)] # filtering
		# matplot(t(var.xy.ms), t='l', lty=1, col=c(1,3))
		var.trajectories[[i.theta]] <- var.xy.ms


		# spiking - from trajectory + ratemaps
		# t0 <- proc.time()[1]
		if (code=='MAP')	spikes <- sim.spikes.x(xy.ms, xy.ms, pcs.rate, rep(1, T.theta), dt=0.001, gain.spw=1)
		if (code=='sampling')	spikes <- sim.spikes.x(xy.ms, xy.ms, pcs.rate, rep(1, T.theta), dt=0.001, gain.spw=1)
		if (code=='PPC') spikes <- sim.spikes.x.ppc(xy.ms, var.xy.ms, pcs.rate, dt=0.001, conv.enc=F)
		if (code=='PPC_var') spikes <- sim.spikes.x.ppc(xy.ms, var.xy.ms, pcs.rate, dt=0.001, conv.enc=F, use_SD=F, sigma0=0.0025)
		if (code=='DDC')	spikes <- sim.spikes.x.ppc(xy.ms, var.xy.ms, pcs.rate, dt=0.001, conv.enc=T)
		# t1 <- proc.time()[1]
		# print (t1-t0)
		
		# collecting spikes into a matrix of spike times and cell ids.
		cell.spikes <- raster2spt(spikes[i_sort,])	
		phi.t <- cell.spikes[,1] / (T.theta/1000) * 2 * pi
		cell.spikes[,1] <- cell.spikes[,1] + t0
		cell.spikes <- cbind(cell.spikes, phi.t)
		spt <- rbind(spt, cell.spikes)
		
		
		t0 <- t0 + T.theta/1000
		theta.starts <- c(theta.starts, t0)
		if (t0 > reftau ) {
			cat(t0, ' ')
			reftau <- reftau + 10
		}
	}
	
	cat('spikes generated... \n')		
	spt <- spt[-1,]
	dimnames(spt) <- list(NULL, c("time", "cell", "phi"))
	
	############################################################
	## writing out data
	############################################################

	data <- list(spt=spt, # spikes
		theta.starts=theta.starts,
		xy.trajectories=xy.trajectories,
		var.trajectories=var.trajectories,
		rates=rates) # place cell firing rates
		
 	outfile <- paste('spikes1D_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, '.RData', sep='')
	save(data, file=outfile)
}


############################################################
## theta sequences and phase precession in single trial
############################################################
 code <- 'PPC_var'
 N.cells <- 200
 spPth <- 40
 Tmax <- 396.7
 
motionfile <- paste('motiondata1D_Tmax', Tmax, '.RData', sep='')
load(motionfile) 
spikefile <- paste('spikes1D_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, '.RData', sep='')
load(spikefile)

T.axis <- motiondata$T.axis
xy <- motiondata$position
spt <- data$spt


plot(T.axis, xy$x[1,], t='l', xlim=c(298, 302))
T_min <- 300.5
T_max <- 301
tt <- seq(T_min, T_max, by=1/1000)

TseqFile <- paste('SimFigs/09_1Dsim/Tseqs_', code, '_T', T_min, '.pdf', sep='')
pdf(file=TseqFile, 5, 3, useD=F)
jj <- which((spt[,1] > T_min) & (spt[,1] < T_max))
plot(spt[jj,1],spt[jj, 2], pch=16, col=viridis(200)[spt[jj, 2]], axes=F, xlab='time (s)', ylab='cells', ylim=c(0, 200))
lines(tt, 50*cos(tt*2*pi*10) + 100)
axis(1, seq(T_min, T_max, by=0.1), round(seq(T_min, T_max, by=0.1) - T_min, 2))
axis(2, las=2)
for (i in 1:5) lines(1:100/1000+T_min+(i-1)*0.1, data$xy.trajectories[[2994+i]][1,]*100, lwd=3, col=grey(0.75))

dev.off()

#############################################################
## across trials
##############################################################
fxpos <- xy$x[1,]
TT <- length(fxpos)

th_right <- 1.98
th_left <- 0.02

istart.right <- which((fxpos[-TT] < th_left) & (fxpos[-1] > th_left))
iend.right <- which((fxpos[-TT] < th_right) & (fxpos[-1] > th_right))
istart.left <- which((fxpos[-TT] > th_right) & (fxpos[-1] < th_right))
iend.left <- which((fxpos[-TT] > th_left) & (fxpos[-1] < th_left))

right.runs <- matrix(c(0,0), 1)
ii_run <- 1
for (i.run in 1:length(istart.right)){
	istart <- istart.right[i.run]
	ii <- which(iend.right > istart)
	if (length(ii) > 0){
		iend <- iend.right[min(ii)]
		if ((iend - istart) <= 150){
			N <- nrow(right.runs) # check if it returned to start ...
			if (iend == right.runs[N,2]){
				right.runs[N,] <- c(istart, iend)
			} else {
				right.runs <- rbind(right.runs, c(istart, iend))
			}
		}
	}
}
right.runs <- right.runs[-1,]

left.runs <- matrix(c(0,0), 1)
for (i.run in 1:length(istart.left)){
	istart <- istart.left[i.run]
	ii <- which(iend.left > istart)
	if (length(ii) > 0){
		iend <- iend.left[min(ii)]
		if ((iend - istart) <= 150){
			N <- nrow(left.runs) # check if it returned to start ...
			if (iend == left.runs[N,2]){
				left.runs[N,] <- c(istart, iend)
			} else {
				left.runs <- rbind(left.runs, c(istart, iend))
			}
		}
	}
}
left.runs <- left.runs[-1,]

PhPrecFile <- paste('SimFigs/09_1Dsim/PhasePrec_', code, '.pdf', sep='')
pdf(file=PhPrecFile, 10, 2, useD=F)
par(mfcol=c(2, 5))
par(mar=c(1,1,2,1))
cellids <- c(65, 77, 80, 83, 91, 112, 117, 123, 143, 148)
for (cellid in cellids){
	xphi <- c(0,0)
	# for (irun in 1:nrow(left.runs)){
		# T_min <- T.axis[left.runs[irun,1]]
		# T_max <- T.axis[left.runs[irun,2]]
	for (irun in 1:nrow(right.runs)){
		T_min <- T.axis[right.runs[irun,1]]
		T_max <- T.axis[right.runs[irun,2]]
		
		jj <- which((spt[,1] > T_min) & (spt[,1] < T_max) & (spt[,2] == cellid))
		t_spikes <- spt[jj, 1]
		phi_spikes <- spt[jj, 3]
		x_spikes <- approx(T.axis, xy$x[1,], t_spikes)$y
		xphi <- rbind(xphi, cbind(x_spikes, phi_spikes))
	}
	
	xphi <- xphi[-1,]
	plot(xphi[,1], xphi[,2], xlim=c(0, 2), ylim=c(0, 2*pi), pch=16, col=grey(0, alpha=0.25), cex=1, main=cellid, xlab='', ylab='', axes=F)
	rect(0, 0, 2, 2*pi)
	# readline(cellid)
}
dev.off()

