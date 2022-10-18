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
require(Matrix)
require(viridis)
require(gplots)
require(ellipse)
ell2 <- ellipse::ellipse
source('~/Projects/PATTERN/ThetaSeq/gen_inf_kalman.R', chdir = TRUE)
source('~/Programs/Code/HPC_Data/SimRunND.R', chdir = TRUE)
source('~/Projects/PATTERN/ThetaSeq/FosterData/PlaceCellFunctions.R', chdir = TRUE)
source('~/Projects/PATTERN/ThetaSeq/FosterData/SpeedFunctions.R', chdir = TRUE)

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

Tmax <- 30*60 # 5 min
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
vx <- sim.run(Tmax=Tmax, xmax=L, dt=dt.run, pars.speed, ndim=2, seed=11) 

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
if (sim.irreg) outfile <- paste('motiondata_Tmax', Tmax, 'irreg.RData', sep='') else outfile <- paste('motiondata_Tmax', Tmax, '.RData', sep='')
if (outfile %in% list.files()){
	load(file=outfile)
} else {
	post <- local.infer.kalman(xy$y, q, r, drive, N.samples, pars.speed, X=xy$x, graphics=F, irreg=sim.irreg)
	# each theta cycle represents the inferred trajectory up to now
	# the last theta cycle starts when sensory data ends
	motiondata <- list(T.axis=T.axis, intended=vx, position=xy, post=post, pars.speed=pars.speed, q=q, r=r, drive=drive)
	save(motiondata, file=outfile)
}


tt <- motiondata$T.axis
xpos <- motiondata$position$x[1,] * 100
ypos <- motiondata$position$x[2,] * 100
# motion_dir <- atan2(diff(ypos), diff(xpos))
# motion_dir <- c(motion_dir[1], motion_dir)
# pos <- cbind(tt, xpos, ypos, motion_dir) 
# txyvd <- posfilter(pos, sdfilt=0.25, graphics=T)
# stats <- plot_motion_stats(xyvd= txyvd[,2:5], dt_pos=dt.run)


tt_30Hz <- seq(1/30, Tmax, by=1/30)
xpos30 <- approx(tt, xpos, tt_30Hz, rule=2)$y
ypos30 <- approx(tt, ypos, tt_30Hz, rule=2)$y
motion_dir30 <- atan2(diff(ypos30), diff(xpos30))
motion_dir30 <- c(motion_dir30[1], motion_dir30)
pos30 <- cbind(tt_30Hz, xpos30, ypos30, motion_dir30) 
txyvd30 <- posfilter(pos30, sdfilt=0.25, graphics=F)
# stats <- plot_motion_stats(xyvd= txyvd30[,2:5], dt_pos=1/30, file='model_motion_stats_new.pdf')


####################################################
## How to select SIMILAR or DIFFERENT samples on subsequent cycles?
## 1. We rotate each sample trajectory by the current running direction and subtract the current position
## 2. We select a start sample trajectory at t=0
## 3. For all 100 proposed trajectories we calculate the difference between the direction of the last sampled trajectory and the current proposal
## 4. We order proposals by their similarity
## 5. We select probabilistically over- or underweighting similar proposals
####################################################
# stats <- plot_motion_stats(xyvd= txyvd30[,2:5], dt_pos=dt.run)

taus <- motiondata$post[[1]]$tstat[1] # theta cycle times
for (tt in 2:length(motiondata$post)) taus <- c(taus, motiondata$post[[tt]]$tstat[1])


rotateX <- function(X, phi, origin=c(0,0), offset=F){
	# 2 x N matrix of position coordinates to rotate
	# phi: rotation angle in radians
	# origin: the center of the rotation
	# shift: if true, then the origin is added to the rotated coordinates

	phi <- as.vector(phi)	
	RotM <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2, 2)
	X0 <- rbind(X[1,] - origin[1], X[2,] - origin[2])

	Xrot <- RotM %*% X0
	if (offset) {
		XX <- rbind(Xrot[1,] + origin[1], Xrot[2,] + origin[2])
	} else {
		XX <- Xrot
	}
	XX
}

sigm <- function(x, ampl=1, slope=1, threshold=0){
	fx <- ampl / (1 + exp(-slope*(x-threshold)))
	fx
}


# pdf(file='intended_motion_sens_post.pdf', 6, 6, useD=F)
# for (i in 1:100){
	ii.theta <- seq(203, 212)+18*10
	t.theta <- taus[ii.theta]
	xpos.theta <- approx(txyvd30[,1], txyvd30[,2], t.theta, rule=2)$y
	ypos.theta <- approx(txyvd30[,1], txyvd30[,3], t.theta, rule=2)$y
	motion_dir.theta <- approx(txyvd30[,1], txyvd30[,5], t.theta, rule=2)$y
	# plot(xpos.theta, ypos.theta, t='o', pch=16, axes=T, xlab='', ylab='', xlim=c(0,200), ylim=c(0,200), lwd=1, col=viridis(10), main=maintext)
	# readline(i)
# }
t.samples <- seq(1, 100, length=21)
t.samples.ms <- seq(1, 100)
Nfilt <- 24
filt <- rep(1, Nfilt)/Nfilt

pars <- cbind(slope=c(5, -5, 5, 8, -8, 0), th=c(pi/2, pi/4, pi/4, pi/8, pi/8, pi/2))
rownames(pars) <- list('strong-anti', 'similar', 'anti-half', 'anti-quarter', 'same', 'uniform')
# beta <- 5
# for (beta in c(5, 3, 1, 0, -1, -3)){ # 0: uniform sampling, 5: cycling, -1: correlated samples
i_par <- 1
for (i_par in 1:6){

	theta_trajectories <- list()
	i.theta <- ii.theta[1]
	j <- 1
	i <- 1
	sam <- motiondata$post[[i.theta]]$samples[i,,]
	Rsam <- rotateX(sam, (-1) * motion_dir.theta[j], origin=c(xpos.theta[j]/100, ypos.theta[j]/100), offset=F)
	dir_past <- coord2rad(Rsam[1,21],Rsam[2,21])
	theta_trajectories[[j]] <- sam
	
	outname <- paste('SimFigs/10_CorrelatedSamples/Sampled_trajectories_', rownames(pars)[i_par], '_i_theta', i.theta, '.pdf', sep='')
	pdf(file=outname, 4, 4, useD=F)
	par(mar=c(1,1,1,1))
	par(mfcol=c(1,1))
	maintext = paste('sampling with ', rownames(pars)[i_par], sep='')
	plot(xpos.theta, ypos.theta, t='o', pch=16, axes=F, xlab='', ylab='', xlim=c(120,220), ylim=c(30,130), lwd=1, col=viridis(10), main=maintext)
	# plot(xpos.theta, ypos.theta, t='o', pch=16, axes=T, xlab='', ylab='', xlim=c(0,200), ylim=c(0,200), lwd=1, col=viridis(10), main=maintext)
	
	sam.ms <- rbind(approx(t.samples, sam[1,], t.samples.ms)$y, approx(t.samples, sam[2,], t.samples.ms)$y) # interpolating
	sam.ms <- cbind(matrix(rep(sam.ms[,1], Nfilt/2), 2), sam.ms, matrix(rep(sam.ms[,100], Nfilt/2), 2)) # padding for filtering
	sam.ms <- rbind(filter(sam.ms[1,], filt), filter(sam.ms[2,], filt))[,(Nfilt/2):(100+Nfilt/2-1)] # filtering
	
	lines(sam.ms[1,]*100, sam.ms[2,]*100, col=viridis(10)[j], t='l')
	points(sam.ms[1,100]*100, sam.ms[2,100]*100, col=viridis(10)[j], t='p', pch=15)
	
	
	for (j in 2:10){
		i.theta <- ii.theta[j]
		err_dir <- rep(NA, 100)
		for (i in 1:100){
			sam <- motiondata$post[[i.theta]]$samples[i,,]
			Rsam <- rotateX(sam, (-1) * motion_dir.theta[j], origin=c(xpos.theta[j]/100, ypos.theta[j]/100), offset=F)
			# points(sam[1,21], sam[2,21], col=grey(0.75))
			dir <- coord2rad(Rsam[1,21],Rsam[2,21])
			err_dir[i] <- min((dir - dir_past) %% (2*pi), (dir_past - dir) %% (2*pi))
		}
		# norm_err_dir <- (err_dir - mean(err_dir)) / sd(err_dir)
		# p_traj <- exp(beta * norm_err_dir)
		p_traj <- sigm(err_dir, slope=pars[i_par,1], th=pars[i_par,2])
		# plot(err_dir, p_traj)
		# print(log10(sort(p_traj, decreasing=T)[1:3]))
		i_traj <- sample(1:100, 1, prob = p_traj)
		# print(i_traj)
		
		sam <- motiondata$post[[i.theta]]$samples[i_traj,,]
		Rsam <- rotateX(sam, (-1) * motion_dir.theta[j], origin=c(xpos.theta[j]/100, ypos.theta[j]/100), offset=F)
		dir_past <- coord2rad(Rsam[1,21],Rsam[2,21])
		theta_trajectories[[j]] <- sam
		
		sam.ms <- rbind(approx(t.samples, sam[1,], t.samples.ms)$y, approx(t.samples, sam[2,], t.samples.ms)$y) # interpolating
		sam.ms <- cbind(matrix(rep(sam.ms[,1], Nfilt/2), 2), sam.ms, matrix(rep(sam.ms[,100], Nfilt/2), 2)) # padding for filtering
		sam.ms <- rbind(filter(sam.ms[1,], filt), filter(sam.ms[2,], filt))[,(Nfilt/2):(100+Nfilt/2-1)] # filtering
	
		lines(sam.ms[1,]*100, sam.ms[2,]*100, col=viridis(10)[j], t='l')
		points(sam.ms[1,100]*100, sam.ms[2,100]*100, col=viridis(10)[j], t='p', pch=15)
		# lines(sam.ms[1,], sam.ms[2,], col=2, t='l')
	}


	scalebar2(20, 0, '20 cm', pos='topleft')
	dev.off()
}


###############################################
## Generate synthetic neuronal data to encode theta sequences

## 1. generate place cell tuning
N.cells <- 200
cat("generating the rate matrix... \n")
# png(file='Figs/placecells16.png', 800, 800)
## Mixture of Gaussian place fields: more realistic, but inference is slower - can only be implemented using the histogram method
pcs.rate <- gen.rate(N.cells, 2, 100, 1, seed=19, active.cells=T, rate.max=c(5, 15),  rate0=c(0.1,0.25), r.place.field=c(0.1, 0.3), graphics=F, rand.nfields=T)
## Multivariate Gaussian place fields: inference can be done analytically
# pcs.rate <- gen.rate(N.cells, 2, 81, 1, seed=19, active.cells=T, rate.max=c(5, 30),  rate0=c(0.1,0.25), r.place.field=c(0.1, 0.344), graphics=F, rand.nfields=F)
# dev.off()

rates <- array(NA, dim=c(N.cells, 40, 40))
xx <- rep(seq(2.5, 197.5, by=5)/100, each=40)
yy <- rep(seq(2.5, 197.5, by=5)/100, 40)
for (i in 1:N.cells) rates[i,,] <- matrix(get.rate.x(xx=rbind(xx, yy), pcs.rate[[i]]), 40, byrow=T)
# hist(apply(rates, 1, mean), br=seq(0, 30, by=2))
# hist(apply(rates*0.1, c(2,3), sum))
spPth <- round(mean(apply(rates*0.1, c(2,3), sum))) # spikes per theta - should be 31
cat('spikes per theta: ', spPth, '±', sd(apply(rates*0.1, c(2,3), sum)))

attr(rates, 'xcenters') <- seq(2.5, 197.5, by=5)
attr(rates, 'ycenters') <- seq(2.5, 197.5, by=5)

# save(rates, file='SimFigs/03DDC_estVars/rates_synth_orig.RData')

# plot_rate_stats(rates, 1:N.cells, 1:N.cells, c(9, 10, 18, 19))


## 1. generate place cell ratemaps
t.samples <- seq(1, 100, length=21)
Nfilt <- 24
filt <- rep(1, Nfilt)/Nfilt

i.code <- 2
codes <- c('MAP', 'sampling', 'PPC', 'DDC', 'PPC_var')
# beta <- 5 # controls autocorrelation of samples. beta > 0: cycling; beta = 0: independent; beta < 0: correlated

N.theta <- length(motiondata$post)
# N.theta <- 1000
xpos.theta <- approx(txyvd30[,1], txyvd30[,2], taus, rule=2)$y
ypos.theta <- approx(txyvd30[,1], txyvd30[,3], taus, rule=2)$y
motion_dir.theta <- approx(txyvd30[,1], txyvd30[,5], taus, rule=2)$y

# for (i.code in 2:2){
i.code <- 2
for (i_par in 1:6){

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
			if (i.theta == 1){
				sam <- motiondata$post[[i.theta]]$samples[1,,] # for sampling	
				Rsam <- rotateX(sam, (-1) * motion_dir.theta[i.theta], origin=c(xpos.theta[i.theta]/100, ypos.theta[i.theta]/100), offset=F)
				dir_past <- coord2rad(Rsam[1,21],Rsam[2,21])
			} else {
				err_dir <- rep(NA, 100)
				for (i.sam in 1:100){
					sam <- motiondata$post[[i.theta]]$samples[i.sam,,]
					Rsam <- rotateX(sam, (-1) * motion_dir.theta[i.theta], origin=c(xpos.theta[i.theta]/100, ypos.theta[i.theta]/100), offset=F)
					dir <- coord2rad(Rsam[1,21],Rsam[2,21])
					err_dir[i.sam] <- min((dir - dir_past) %% (2*pi), (dir_past - dir) %% (2*pi))
				}
				# norm_err_dir <- (err_dir - mean(err_dir)) / sd(err_dir)
				# p_traj <- exp(beta * norm_err_dir)
				p_traj <- sigm(err_dir, slope=pars[i_par,1], th=pars[i_par,2])
				# print(log10(sort(p_traj, decreasing=T)[1:3]))
				i_traj <- sample(1:100, 1, prob = p_traj)
				# print(i_traj)

				sam <- motiondata$post[[i.theta]]$samples[i_traj,,]
				Rsam <- rotateX(sam, (-1) * motion_dir.theta[i.theta], origin=c(xpos.theta[i.theta]/100, ypos.theta[i.theta]/100), offset=F)
				dir_past <- coord2rad(Rsam[1,21],Rsam[2,21])
			}
			xy.tau <- sam
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
		cell.spikes <- raster2spt(spikes)	
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
		
	if (sim.irreg) outfile <- paste('spikes_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, 'irreg.RData', sep='') else outfile <- paste('spikes_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, '_', rownames(pars)[i_par], '.RData', sep='')
	save(data, file=outfile)
}

