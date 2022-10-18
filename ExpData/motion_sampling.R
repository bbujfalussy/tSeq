########################################################################
## Script to reproduce Figure 7–Figure supplement 2. Model-free replication of the main findings of the paper. 
## Generating potential trajectories and drive place cell activity accordingly

library(R.matlab)
library(colormap)
library(viridis)
library(fmsb)
library(ellipse)
draw_ellipse <- ellipse::ellipse

sampling.rate <- 1250 # Hz

tryCatch({data <- readMat("Pfeiffer13/spikes.mat")}, error = function(e) print('File containing the spikes not found. Please contact the authors to obtain it!'), warning=function(w) print('File containing the spikes not found. Please contact the authors to obtain it!'))
data <- data[[1]]
names(data) <- list('rat1', 'rat2', 'rat3', 'rat4')
summary(data)
for (i in 1:4){
	names(data[[i]]) <- list('day1', 'day2')
	for (j in 1:2) names(data[[i]][[j]]) <- list('spikes', 'position', 'I_Cells', 'E_Cells')
	
} 

######################################################
source('./SpeedFunctions.R')

ddata <- data$rat1$day2
ncells <- max(ddata$spikes[,2])

i_rat <- 1
i_day <- 2
cat('loading data for rat', i_rat, 'day', i_day, '... \n')

ddata <- data[[i_rat]][[i_day]]
ncells <- max(ddata$spikes[,2])


## filter position to remove noise and make it smoother
txyvd <- posfilter(ddata$position, sdfilt=0.25, graphics=T)

## detect periods of movement - at least 1s long move with at least 0.2 s stops between
dt_pos <- median(diff(txyvd[,1]))
vmin <- 5 # m/s
T_run_min <- 1 # s
T_stop_min <- 0.2 # s
ii_run <- txyvd[,4] > vmin
ii_run[is.na(ii_run)] <- F
i_run <- get_long_periods(ii_run, T_run_min / dt_pos, T_stop_min / dt_pos)

txyvd[!i_run,4:5] <- NA
		

# ###########################################################
stats <- plot_motion_stats(xyvd= txyvd[,2:5], dt_pos= dt_pos)
		
#########################################################
## speed-conditional trajectories 

N_traj <- 10000
T_future <- 1.5 # second
N_future <- round(T_future / dt_pos)
traj <- array(NA, dim=c(N_traj, 10, N_future+1, 2))
speed_q <- quantile(txyvd[,4], c(0:10/10), na.rm=T)
speed_q[11] <- speed_q[11] + 1
ns <- rep(0, 10)
for (i.time in 1:(nrow(txyvd)-(N_future + 1))){
	traj.i <- txyvd[i.time:(i.time+ N_future),]
	if (sum(is.na(traj.i[,4])) == 0){
		start_speed <- traj.i[1,4]
		i_speed <- max(which(start_speed >= speed_q))
		if (ns[i_speed] < N_traj){
			
			xy <- t(traj.i[,2:3])
			xy <- xy - xy[,1]
			theta <- -traj.i[1,5]
			rotation_matrix <- matrix(c(cos(theta), -sin(theta)	, sin(theta), cos(theta)), 2, 2, byrow=T)
			# rotation_matrix <- matrix(c(1,0,0,1), 2, 2, byrow=T)
			xy_rot <- rotation_matrix %*% xy
			# plot(xy[1,], xy[2,], xlim=c(-10, 20), ylim=c(-15, 15))
			# points(xy_rot[1,], xy_rot[2,], col=2)
			ns[i_speed] <- ns[i_speed] + 1
			traj[ns[i_speed], i_speed, , ] <- t(xy_rot)
		}		
		
	}
}
				
		
par(mfcol=c(2,2))
		
for (i_speed in c(1, 3, 5, 10)){
	ii <- seq(1, min(ns[i_speed], 2000), by=20)
	
	plot(t(traj[i,i_speed,,1]), t(traj[i,i_speed,,2]), t='l', col=grey(0.8), xlim=c(-20, 60), ylim=c(-30, 30), xlab='x coordinate (cm)', ylab='y coordinate (cm)', main=paste('speed quantile', i_speed), axes=F)
	for (i in ii[-1]){
		lines(t(traj[i,i_speed,,1]), t(traj[i,i_speed,,2]), t='l', col=grey(0.8))
	}
	lines(draw_ellipse(diag(c(1, 1)), scale=c(sd(traj[ii,i_speed,N_future+1,1]), sd(traj[ii,i_speed,N_future+1,2])), centre=colMeans(traj[ii,i_speed,N_future+1,]), level=0.75), col=2, t='l')
	axis(1)
	axis(2, las=2)				
}

############################################################################################
### for a few particular position and speed, collect 100 potential future trajectories and compare it with the actual trajectory
###
### 100 potential future trajectory:
###		- matching the current position and speed
### 		- remaining within the boundaries of the box
t0 <- 5000
dir.create('./DataFigures/Fig7_FS2/', showW=F)
png(filename=paste('./DataFigures/Fig7_FS2/sampled_trajectories_', t0, '_1500ms.png', sep=''), 1200, 1200, pointsize='24')
# pdf(file=paste('sampled_trajectories_', t0, '_1500ms.pdf', sep=''), 6, 6, useD=F)
par(mar=c(1,1,1,1))
n_plotted <- 0
for (i.time in (t0+seq(1, 2005, by=100))){
	traj.i <- txyvd[i.time:(i.time+ N_future),]
	if (sum(is.na(traj.i[,4])) == 0){
		start_speed <- traj.i[1,4]
		i_speed <- max(which(start_speed >= speed_q))
			
		dist_traj <- traj[,i_speed,,]
		xy <- t(traj.i[,2:3])
		theta <- traj.i[1,5]
		rotation_matrix <- matrix(c(cos(theta), -sin(theta)	, sin(theta), cos(theta)), 2, 2, byrow=T)
		accepted_traj <- array(NA, dim=c(100, N_future+1, 2))
		n_accepted <- 0
		
		while (n_accepted < 100){
			i <- sample(1:ns[i_speed], 1)
			candidate <- dist_traj[i,, ]
			candidate_rot <- rotation_matrix %*% t(candidate) +  xy[,2]
			if ((max(candidate_rot) < 200)&(min(candidate_rot) > 0)){
				n_accepted <- n_accepted + 1
				accepted_traj[n_accepted,,] <- t(candidate_rot)
			}
		}

		n_plotted <- n_plotted + 1
		if (n_plotted == 1){	
			plot(accepted_traj[1,,1], accepted_traj[1,,2], t='l', col=grey(0.75), xlim=c(0, 200), ylim=c(0,200), axes=F, xlab='', ylab='')
		}
		for (i in 2:100){
			lines(accepted_traj[i,,1], accepted_traj[i,,2], t='l', col=grey(0.75), xlim=c(0, 200), ylim=c(0,200), axes=F, xlab='x', ylab='y')
		}
		lines(xy[1,], xy[2,], t='l')		
	}
}

rect(0, 0, 200, 200, lwd=3)
dev.off()

####################################################################################
## collect data to drive theta sequences
T.axis <- NA
position <- list(x=c(NA, NA, NA, NA))
post <- list()
dt_orig <- median(diff(txyvd[,1]), na.rm=T)
dt_new <- 0.05
T0 <- min(txyvd[,1], na.rm=T)
tsam_orig <- seq(0, by=dt_orig, length=N_future+1)
tsam_new <- seq(0, by=dt_new, length=31)
k <- 1

for (i.time in (seq(1,nrow(txyvd)-N_future, by=3))){
	traj.i <- txyvd[i.time:(i.time+ N_future),]
	T.axis <- c(T.axis, traj.i[1,1] - T0)
	position$x <- rbind(position$x, traj.i[1,2:5])
	if (sum(is.na(traj.i[,4])) == 0){
		start_speed <- traj.i[1,4]
		theta <- traj.i[1,5] # the direction and the speed must be calculated using the previous two positions
		
		i_speed <- max(which(start_speed >= speed_q))
		dist_traj <- traj[,i_speed,,]
		xy <- t(traj.i[,2:3])
		rotation_matrix <- matrix(c(cos(theta), -sin(theta)	, sin(theta), cos(theta)), 2, 2, byrow=T)
		samples.tau <- array(NA, dim=c(100, 2, 31))
		n_accepted <- 0
		
		while (n_accepted < 100){
			i <- sample(1:ns[i_speed], 1)
			candidate <- dist_traj[i,, ]
			candidate_rot <- rotation_matrix %*% t(candidate) +  xy[,2]
			if ((max(candidate_rot) < 200)&(min(candidate_rot) > 0)){
				n_accepted <- n_accepted + 1
				x_sam <- approx(tsam_orig, candidate_rot[1,], tsam_new)$y
				y_sam <- approx(tsam_orig, candidate_rot[2,], tsam_new)$y
				
				samples.tau[n_accepted,,] <- rbind(x_sam, y_sam)
			}
		}

		estimpars.tau <- list(mu=t(apply(samples.tau, c(2,3), mean)), var=t(apply(samples.tau, c(2,3), sd)^2))
		post.tau <- list()
		post.tau$samples <- samples.tau
		post.tau$estimpars <- estimpars.tau
		post.tau$tstats <- list(tau=traj.i[1,1]-T0, t.start=0, t.end=1, t.theta=0.1)
		post[[k]] <- post.tau	
	}
	k <- k + 1
	if (k %% 100==0) print(k)
}

T.axis <- T.axis[-1]
position$x <- position$x[-1,]

sampled_motion <- list(T.axis=T.axis, position=position, post=post)
summary(sampled_motion)
filename <- 'Rat1Day2_t1500_sampled_motion_post.RData'
save(sampled_motion, file=filename)



###########################################################################################
#### generating neuronal activity - synthetic neuronal data to encode theta sequences
###########################################################################################
source('./PlaceCellFunctions.R', chdir = TRUE)

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

plot_rate_stats(rates, 1:N.cells, 1:N.cells, c(9, 10, 18, 19))


## 1. generate place cell ratemaps
t.samples <- seq(1, 100, length=21)
Nfilt <- 24
filt <- rep(1, Nfilt)/Nfilt
dt.theta <- median(diff(motiondata$T.axis))

codes <- c('MAP', 'sampling', 'DDC', 'PPC_var')


motiondata <- sampled_motion
N.theta <- length(motiondata$post)
# N.theta <- 1000

for (i.code in c(1,2,3,4)){
	code <- codes[i.code]
	# code <- 'DDC' # sampling, PPC, DDC, MAP
	reftau <- 30
	spt <- c(0, 0, 0, 0) # spike times, cellids, spike phases, theta cycle
	xy.trajectories <- list() # the smooth trajectory, used for generating place cell activity
	var.trajectories <- list()
	
	t0 <- 0 # start of the current theta cycle
	theta.starts <- t0 # vector with the theta start times

	for (i.theta in 1:N.theta) { # 191: time after which the smoothing achieved steady state?
		if (length(motiondata$post[[i.theta]]) > 0){
			xy.tau <- t(motiondata$post[[i.theta]]$estimpar$mu) / 100 # for MAP, DDCal and PPC, the posterior in "motion-coordinates"
			var.xy.tau <- t(motiondata$post[[i.theta]]$estimpar$var) / 100 / 100 + 0.0016 # for PPC - we add 4cm SD, otherwise PPC would diverge...
			if (code=='sampling')	{
				xy.tau <- motiondata$post[[i.theta]]$samples[1,,] / 100 # for sampling	
			} 
			L.theta <- ncol(xy.tau)
			
			T.theta <- motiondata$post[[i.theta]]$tstats$t.theta * 1000 # this is the real duration of the theta sequence in ms
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
			if (code=='PPC_var') spikes <- sim.spikes.x.ppc(xy.ms, var.xy.ms, pcs.rate, dt=0.001, conv.enc=F, use_SD=F, sigma0=0.0025)
			if (code=='DDC')	spikes <- sim.spikes.x.ppc(xy.ms, var.xy.ms, pcs.rate, dt=0.001, conv.enc=T)
			# t1 <- proc.time()[1]
			# print (t1-t0)
			
			# collecting spikes into a matrix of spike times and cell ids.
			cell.spikes <- raster2spt(spikes)	
			phi.t <- cell.spikes[,1] / (T.theta/1000) * 2 * pi
			spt_theta_cycle <- rep(i.theta, nrow(cell.spikes))
			cell.spikes[,1] <- cell.spikes[,1] + t0
			cell.spikes <- cbind(cell.spikes, phi.t, spt_theta_cycle)
			spt <- rbind(spt, cell.spikes)
		} else {
			T.theta <- dt.theta
		}
		
		t0 <- motiondata$T.axis[i.theta]  # start of the current theta cycle
		theta.starts <- c(theta.starts, t0) # vector with the theta start times
		if (t0 > reftau ) {
			cat(t0, ' ')
			reftau <- reftau + 10
		}
	}

	cat('spikes generated... \n')		
	spt <- spt[-1,]
	theta.starts <- theta.starts[-1]
	dimnames(spt) <- list(NULL, c("time", "cell", "phi", 'cycle'))
		
	############################################################
	## writing out data
	############################################################

	data <- list(spt=spt, # spikes
		theta.starts=theta.starts,
		xy.trajectories=xy.trajectories,
		var.trajectories=var.trajectories,
		rates=rates) # place cell firing rates
	
	Tmax <- t0
	outfile <- paste('../Fig7_SF2/spikes_Tmax', round(Tmax), '_Ncells', N.cells, '_', code, '_spPth', spPth, '.RData', sep='')
	save(data, file=outfile)
}


