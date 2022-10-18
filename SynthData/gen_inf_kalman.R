############################################
##############################################
## functions to generate data and perform inference under the Kalman filter model 
## detailed in Figure 3–Figure supplement 1a. 
## In short, the agent is trying to follow an extarnally set goal path, 
## using noisy motor commands and its noisy position estimate
## x is the position
## y is the observation
## r is the variance of the observation noise
## q is the variance of the motor noise
## u is the motor command (drive)


source('./SimRunND.R', chdir = TRUE)
source('./utils.R', chdir = TRUE)
require(abind)
require(ellipse)

sim.kalman.step <- function(x, q=1, r=1, u=0){
## simulating a single Kalman step
## q is the evolution noise variance
## r is the observation noise variance
## u is the drive
	n <- length(x)
	xx <- x + u + rnorm(n, 0, sqrt(q))
	yy <- xx + rnorm(n, 0, sqrt(r))
	list(x=xx, y=yy)	
}


##############################################
## 3. Perform the inference
## 3.1 the forward probabilities
## see 16.3 in MurphyMLPP.pdf
filt.kalman <- function(obs, q, r, prior, U=NULL){
	# calculates the probability = P(X_t | Y_1:t)
	# providing the prior, P(X_t-1 | Y_1:t-1) it can be used to update the posterior by a single step
	# inference has a particularly simple form, since we assumed that the noise is diagonal and A and C are identity
	# the posterior variance is therefore diagonal
	# NA in the observation means missing observation, causing prediction only!
	
	n <- nrow(obs) # number of variables
	Tmax <- ncol(obs) # number of time points
	post.mu <- matrix(NA, n, Tmax)
	post.var <- matrix(NA, n, Tmax)

	if (is.null(U)) U <- matrix(0, n, Tmax)
	if (is.null(prior))	prior <- list(mu=obs[,1], var=c(r, r)	)
	
	prior.mu <- prior$mu
	prior.var <- prior$var

	for (t in 1:Tmax){
		mu.pred <- prior.mu + U[,t]
		var.pred <- prior.var + q

		is.obs <- prod(!is.na(obs[,t]))
		if (is.obs) { # observation is present
			e.t <- obs[,t] - mu.pred
			K <- var.pred / (var.pred + r)
	
			post.mu[,t] <- mu.pred + K * e.t
			post.var[,t] <- (1-K) * var.pred
		
			prior.mu <- post.mu[,t]
			prior.var <- post.var[,t]
		} else {
			post.mu[,t] <- mu.pred
			post.var[,t] <- var.pred
			
		}
	}

	post <- list(mu=post.mu, var=post.var)
	post
}


## 3.2 the backward smoothing - based on 16.3.2.1 of MurphyMLPP.pdf
smooth.kalman <- function(filt.post, q, U=NULL){
	# calculates the probability = log(P(Y_t+1:T | X_t))
	fmu <- filt.post$mu
	fvar <- filt.post$var
	if (is.null(U)) U <- matrix(0, n, Tmax)

	n <- nrow(fmu)
	Tmax <- ncol(fmu)

	mu <- matrix(NA, n, Tmax)
	var <- matrix(NA, n, Tmax)

	mu[,Tmax] <- fmu[,Tmax]
	var[,Tmax] <- fvar[,Tmax]

	for (t in rev(1:(Tmax-1))){
		fvar.t <- fvar[,t]
		fvar.tp <- fvar.t + q
		J <- fvar.t / fvar.tp
		
		fmu.tp <- fmu[,t] + U[,t+1]
		mu.t <- fmu[,t] + J * (mu[,t+1] - fmu.tp)
		var.t <- fvar.t + J^2 * (var[,t+1] - fvar.tp)
		
		mu[,t] <- mu.t
		var[,t] <- var.t
	}
	post <- list(mu=mu, var=var)
	post
}


#########################################################
## sampling trajectories from the posterior
## sampling is started at the end of the trajectory proceeding backwards:
## P(x_T-1|X_T, Y_1:T) = P(x_T-1|X_T, Y_1:T-1) = P(x_T-1| Y_1:T-1) * P(X_T|X_T-1)/P(X_T|Y_1:T-1)
## where the first term is given by the filtering distribution
## prediction is returned if required

# n.sam <- 100
# filt.post <- list(mu=filtered.post$mu[,2996:3000], var=filtered.post$var[,2996:3000])
# U <- drive[,2996:3000]
# n.predict <- 50
# dt <- dt.run
# xmax=1

sample.kalman <- function(n.sam, filt.post, q, U=NULL, n.predict=0, pars.speed=NULL){
	## alpha: smoothness factor for calculating the drive
	## planned.pos: position planned for the future - should be a matrix with n rows and n.predict columns
	
	n <- nrow(filt.post$mu) # number of dimensions
	Tmax <- ncol(filt.post$mu)
	if (is.null(U)) U <- matrix(0, n, Tmax)
	sams <- array(NA, dim=c(n.sam, n, Tmax))
	mu.t <- filt.post$mu[,Tmax]
	var.t <- filt.post$var[,Tmax]
	sams[, 1, Tmax] <- rnorm(n.sam, mu.t[1], sqrt(var.t[1]))
	sams[, 2, Tmax] <- rnorm(n.sam, mu.t[2], sqrt(var.t[2]))
	
	if (Tmax > 1){
		for (i.sam in 1:n.sam){
			for (t in rev(1:(Tmax-1))){
				x.tp <- sams[i.sam, ,t+1]
				drive.tp <- U[,t+1]
				
				mu.filt.t <- filt.post$mu[,t]
				var.filt.t <- filt.post$var[,t]
				
				var.t <- 1 / (1/var.filt.t + 1/q)
				mu.t <- var.t * ((x.tp - drive.tp) / q + mu.filt.t/var.filt.t)
	
				sams[i.sam, , t] <- rnorm(2, mu.t, sqrt(var.t))
			}
		}
	}

	if (n.predict > 0){
		sams.pred <- array(NA, dim=c(n.sam, n, n.predict+1))
		sams.pred[, , 1] <- sams[,,Tmax]
		if (is.null(pars.speed)) {
			warning('movement parameters are missing')
			add.drive <- F
		} else {
			add.drive <- T
		}
			
		for (i.sam in 1:n.sam){
			if (add.drive) drive.t <- U[,Tmax]
			velocity.t <- (sams[i.sam,,Tmax] - sams[i.sam,,Tmax-1]) / pars.speed$dt
			speed.t <- sqrt(sum(velocity.t^2))

			for (t in 2:(n.predict+1)){
				last.pos <- sams.pred[i.sam, , t-1]
				if (add.drive) {
					sv <- iterate.sv(speed.t, velocity.t, pars.speed$dt, pars.speed)
					speed.t <- sv$speed
					velocity.t <- sv$velocity
					next.pos <- last.pos + velocity.t * pars.speed$dt
					vx <- check.boundary.weak(velocity.t, next.pos, pars.speed$xmax)
					velocity.t <- vx$v # this does not change the speed only the direction!
					next.pos <- vx$x

					drive.t <- drive.t *	(1-pars.speed$alpha) + pars.speed$alpha * (next.pos - last.pos)
					mu.pred.t <- last.pos + drive.t
				} else {
					mu.pred.t <- last.pos					
				}
				var.pred.t <- q
				sams.pred[i.sam, , t] <- rnorm(2, mu.pred.t, sqrt(var.pred.t))				
			}
		}
		
		sams <- abind(sams, sams.pred[,,-1])
	}

	# matplot(t(sams[,1,]), t='l', lty=1, col=grey(0.6))
	# matplot(t(sams.pred[,1,]), t='l', lty=1, col=grey(0.6))
	# matplot(t(sams[,2,]), t='l', lty=1, col=grey(0.6))
	# matplot(t(sams.pred[,2,]), t='l', lty=1, col=grey(0.6))
	# plot(apply(sams[,2,], 2, sd), xlim=c(0, 16), ylim=c(0, 0.1))
	# points(apply(sams[,1,], 2, sd), col=2)
	# lines(sqrt(smooth.post$var[1,2996:3000]), col=2)

	# points(5:15, apply(sams.pred[,2,], 2, sd), col=3, cex=0.7, pch=16)
	# points(5:15, apply(sams.pred[,1,], 2, sd), col=4, cex=0.7, pch=16)

	# plot(sams[1,1,], sams[1,2,], t='o', pch=16, cex=0.5, xlim=c(-0.1,1.1), ylim=c(-0.1, 1.1), col=viridis(55))
	# for (i in 2:100) points(sams[i,1,], sams[i,2,], t='o', pch=16, cex=0.5, col=viridis(55))
	
	sams
}


## a function that performs online inference and sampling given the observations
## if irreg  == T random past and future intervals and random theta durations 
local.infer.kalman <- function(Y, q, r, drive, N.samples, pars.speed, X=NULL, graphics=F, plot.samples=F, irreg=T){
	# this function calculates the online posterior - P(x_(tau-t.past):tau | y_(1:tau))
	# and samples from the distribution P(x_(tau-t.past):(tau+t.future) | y_(1:tau))
	# Y: observation matrix 2 x L
	# q: evolution noise variance
	# r: observation noise variance
	# drive: evolution drive
	# N.samples: numer of samples to be drawn from the posterior
	# pars.speed: parameters for prediction:
	# X: true position matrix  2 x L
	# irreg: whether each theta cycle should be different or not
	# t.past, t.future: intervals for past and future trajectory duration

	# 	controlling the smoothness of the drive, how the velocity changes in time, the timestep and the borders of the the arena
	# Y <- xy$y; X <- xy$x; graphics <- T; plot.samples <- T
	# IMPORTANT: when irreg == True each cycle is different:
	# - they start at different point in the past -1.5 < t.past < 0 
	# - they move to a different distance into the future 0 < t.future < 1.5
	# - with encoding 1 < t.total < 3 s 
	# - each cycle lasts for 80 < t.theta < 160 ms - so the encoding speed is also variable - built in the local.infer.kalman.irreg function

	dt.run <- 0.1
	Tmax <- ncol(Y) * dt.run # s
	if (irreg == T){
		t.past <- c(-1.5, -0.1); t.future <- c(0, 1.5); t.total <- c(1,3)
		T.past <- t.past / dt.run
		T.future <- t.future / dt.run
		T.total <- t.total / dt.run
		PP <- 11-abs(10:30 - 20) # weights for more frequent selection of typical trajectories
	} else {
		t.past <- c(-1, -1); t.total <- c(2, 2)
		T.past <- t.past / dt.run
	}
	
	if (Tmax < max(t.total)) stop('Data too short, local inference can not be done. Implement full inference instead!')
	if (graphics & !is.null(X)) plot.pos <- T else plot.pos <- F

	post <- list()
	##############################################################
	### first step, initial inference is made	
	tau <- abs(t.past[1]) + dt.run # tau is the current time in ms
	L.filt <- abs(T.past[1]) + 1 # length of the inference part to update in each cycle 
	## we keep the current estimate + 10 from the past and 10 predictions 
	iY <- round(tau / dt.run)# + 1 # position index - same as L.filt for now
	k <- 1 # theta cycle index
	
	
	cat("calculating local inference around tau =", tau, '\n')
	Y.tau <- Y[,1:iY]
	drive.tau <- drive[,1:iY]
	if (!is.null(X)) {
		X.tau <- X[,1:iY]
	} else X.tau <- NULL
	
	prior <- list(mu=Y[,1], var=c(r,r))
	filt.tau <- filt.kalman(Y.tau, q, r, prior, drive.tau)
	prior <- list(mu=filt.tau$mu[, L.filt], var=filt.tau$var[, L.filt]) # the last element of the filtering posterior is going to be the prior for the next step
	smooth.tau <- smooth.kalman(filt.tau, q, U=drive.tau)

	# we take random intervals back and forward. The present is always included, and the duration varies
	if (irreg == T) {
		L.seq <- sample(T.total[1]:T.total[2], 1, prob=PP) # first select the length of the represented trajectory - in steps
		T.start <- sample(T.past[1]:(T.future[2]-L.seq), 1) # we take this many steps back
		if (T.start > -1) T.start <- sample(min(-1, -T.start):(-1), 1) # we need to start from the past - at least 2 steps, to initiate predictions
		T.end <- T.start + L.seq # end = start + length
		if (T.end < 0) T.end <- sample(0:(-T.end), 1) # we want to include some prediction as well
	} else {
		T.start <- -10 # we take this many steps back
		T.end <- 10 # end = start + length
	}
	# N <- 1000
	# Ts <- matrix(NA, 2, N)
	# for (i in 1:N){
		# L.seq <- sample(T.total[1]:T.total[2], 1, prob=PP) # first select the length of the represented trajectory - in steps
		# T.start <- sample(T.past[1]:(T.future[2]-L.seq), 1) # we take this many steps back
		# if (T.start > -1) T.start <- sample(min(-1, -T.start):(-1), 1) # we need to start from the past - at least 2 steps, to initiate predictions
		# T.end <- T.start + L.seq # end = start + length
		# if (T.end < 0) T.end <- sample(0:(-T.end), 1) # we want to include some prediction as well
		# Ts[,i] <- c(T.start, T.end)	
	# }

	# par(mfcol=c(1,3))
	# hist(Ts[1,], br=-16:0+0.5)
	# hist(Ts[2,], br=0:16-0.5)
	# hist(Ts[2,] - Ts[1,], br=10:31-0.5)	
	
	filt.Tau <- list(mu=matrix(filt.tau$mu[,(L.filt + T.start):L.filt], 2), var=matrix(filt.tau$var[,(L.filt + T.start):L.filt], 2))
	smooth.Tau <- list(mu=matrix(smooth.tau$mu[,(L.filt + T.start):L.filt], 2), var=matrix(smooth.tau$var[,(L.filt + T.start):L.filt], 2))
	drive.Tau <- matrix(drive.tau[,(L.filt + T.start):L.filt], 2)
	samples.tau <- sample.kalman(N.samples, filt.Tau, q, U=drive.Tau, n.predict=T.end, pars.speed=pars.speed)
	estimpars.tau <- list(mu=apply(samples.tau, c(2,3), mean), var=apply(samples.tau, c(2,3), var))
	
	if (irreg == T) {
		t.theta.tau <- round(runif(1, 2 * L.seq + 60, 2 * L.seq + 100) / 1000, 3) # s; the duration of the theta cycle encoding this sequence - we will move on this much
	} else {
		t.theta.tau <- 0.1
	}
	post.tau <- list()
	post.tau$smooth <- smooth.Tau
	post.tau$samples <- samples.tau
	post.tau$estimpars <- estimpars.tau
	post.tau$tstats <- c(tau=tau, t.start=T.start*dt.run, t.end=T.end*dt.run, t.theta=t.theta.tau)
	post[[k]] <- post.tau	
	
	if (graphics){
		plot.kalman.samples(post.tau, pars.speed, X.tau[,(L.filt + T.start):L.filt], plot.samples=plot.samples)
	}

	
	#####################################################
	# move on to the next theta cycle
	tau <- tau + t.theta.tau 
	iY.new <- round(tau / dt.run)
	k <- k + 1
	
	while (tau < Tmax){
		cat("calculating local inference around tau =", round(tau, 3), '\n')

		dY <- iY.new - iY # add new sensory input if there is any
		if (dY > 0){
			for (ii in 1:dY) { # we doo the filtering step by step
				Y.ttau <- matrix(Y[,iY+ii], 2) # ttau: conditioned on tau; tau: conditionned on :tau
				drive.ttau <- matrix(drive[,iY+ii], 2) # for filtering

				filt.ttau <- filt.kalman(Y.ttau, q, r, prior, drive.ttau)
				filt.tau$mu <- cbind(filt.tau$mu[,-1], filt.ttau$mu)
				filt.tau$var <- cbind(filt.tau$var[,-1], filt.ttau$var)
			}
			prior <- list(mu=filt.tau$mu[, L.filt], var=filt.tau$var[, L.filt]) # the last element of the filtering posterior is going to be the prior for the next step
		
			drive.tau <- drive[,(iY.new-L.filt+1):iY.new]
			if (!is.null(X)) {
				X.tau <- X[,(iY.new-L.filt+1):iY.new]
			} else X.tau <- NULL	
			smooth.tau <- smooth.kalman(filt.tau, q, U=drive.tau)
		}

		# we take random interwals back and forward. The present is always included, and the duration varies
		if (irreg == T) {
			L.seq <- sample(T.total[1]:T.total[2], 1, prob=PP) # first select the length of the represented trajectory - in steps
			T.start <- sample(T.past[1]:(T.future[2]-L.seq), 1) # we take this many steps back
			if (T.start > -1) T.start <- sample(min(-1, -T.start):(-1), 1) # we need to start from the past - at least 2 steps, to initiate predictions
			T.end <- T.start + L.seq # end = start + length
			if (T.end < 0) T.end <- sample(0:(-T.end), 1) # we want to include some prediction as well
		} else {
			T.start <- -10 # we take this many steps back
			T.end <- 10 # end = start + length			
		}
		
		filt.Tau <- list(mu=matrix(filt.tau$mu[,(L.filt + T.start):L.filt], 2), var=matrix(filt.tau$var[,(L.filt + T.start):L.filt], 2))
		smooth.Tau <- list(mu=matrix(smooth.tau$mu[,(L.filt + T.start):L.filt], 2), var=matrix(smooth.tau$var[,(L.filt + T.start):L.filt], 2))
		drive.Tau <- matrix(drive.tau[,(L.filt + T.start):L.filt], 2)
		samples.tau <- sample.kalman(N.samples, filt.Tau, q, U=drive.Tau, n.predict=T.end, pars.speed=pars.speed)
		estimpars.tau <- list(mu=apply(samples.tau, c(2,3), mean), var=apply(samples.tau, c(2,3), var))

		if (irreg == T) {
			t.theta.tau <- round(runif(1, 2 * L.seq + 60, 2 * L.seq + 100) / 1000, 3) # ms; the duration of the theta cycle encoding this sequence - we will move on this much
		} else {
			t.theta.tau <- 0.1
		}
		post.tau <- list()
		post.tau$smooth <- smooth.Tau
		post.tau$samples <- samples.tau
		post.tau$estimpars <- estimpars.tau
		post.tau$tstats <- c(tau=tau, t.start=T.start*dt.run, t.end=T.end*dt.run, t.theta=t.theta.tau)
		post[[k]] <- post.tau	
		
		if (graphics){
			plot.kalman.samples(post.tau, pars.speed, X.tau[,(L.filt + T.start):L.filt], plot.samples=plot.samples)
		}
		#####################################################
		# move on to the next theta cycle
		iY <- iY.new
		tau <- tau + t.theta.tau 
		iY.new <- round(tau / dt.run)
		k <- k + 1
	}	
	post
}


plot.kalman.samples <- function(post, pars.speed, X=NULL, plot.samples=T){
	smooth <- post$smooth
	samples <- post$samples
	estimpars <- post$estimpars
	
	N.samples <- dim(samples)[1]
	TT <- dim(samples)[3]
	xyrange <- c(-0.2, pars.speed$xmax + 0.2)
	T.past <- ncol(smooth$mu)
	
	
	cols.bg <- viridis(TT, option='B')
	cols.bg.sam <- viridis(TT, option='B', alpha=0.35)
	cols.bg.true <- viridis(TT, option='B')
	# cols.fg <- rgb(1, 0.5, 0)
	# cols.fg.true <- rgb(0, 0.75, 1)
	cols.fg.true <- rgb(128/255, 1, 0/255)
	cols.fg <- viridis(5)[4]

	plot(estimpars$mu[1,], estimpars$mu[2,], xlim=xyrange, ylim=xyrange, t='n', xlab='', ylab='', axes=F)
	scalebar2(0.2, 0, '20 cm', '', 'topleft')
	if (plot.samples){
		for (i in 1:N.samples){
			lines(samples[i,1,], samples[i,2,], t='l', col=grey(0.75, alpha=0.5))
			lines(samples[i,1,], samples[i,2,], t='p', cex=0.8, pch=16, col=cols.bg.sam)
		}
	} else {
		for (i in 1:TT){
			contour <- t(t(ellipse::ellipse(diag(estimpars$var[,i]))) + estimpars$mu[,i])
			lines(contour, col=cols.bg[i])
		}
	}
	# points(smooth$mu[1,], smooth$mu[2,], t='o', cex=1, pch=21, col=cols.fg, bg=cols.bg)
	points(estimpars$mu[1,], estimpars$mu[2,], t='o', cex=1.3, pch=21, col=cols.fg, bg=cols.bg.true, lwd=1)
	points(estimpars$mu[1,T.past], estimpars$mu[2,T.past], cex=1.5, pch=21, col=cols.fg, bg=1, lwd=2)
	if (!is.null(X)) {
		# lines(X[1,], X[2,], t='o', cex=1.3, pch=24, col=cols.fg.true, bg=cols.bg.true, lwd=2)
		lines(X[1,], X[2,], t='l', lwd=2, col=cols.fg.true)
		points(X[1,T.past], X[2,T.past], cex=1.8, pch=24, col= cols.fg.true, bg=1, lwd=2)
	}
	
	# if (plot.samples){
		# legend('topright', legend=c('current', 'samples', 'posterior mean', 'true trajectory'), pt.cex=c(1.8,0.7,1.3,1.3), col=c(cols.fg, cols.bg.sam[round(TT/2)], cols.fg, cols.fg.true), pch=c(21,21,21,24), pt.bg=c(1, cols.bg.sam[round(TT/2)], rep(viridis(3)[2],2)), lwd=c(1,1,2,2))
	# } else {
		# legend('topright', legend=c('present', 'posterior variance', 'posterior mean', 'true trajectory'), pt.cex=c(1.8,1.8,1.3,1.3), col=c(cols.fg, cols.bg[round(TT/2)], cols.fg, cols.fg.true), pch=c(21,21,21,24), pt.bg=c(1, grey(1), rep(viridis(3)[2],2)), lwd=c(2,1,2,2))	
	# }
}


plot.kalman.trajectory <- function(post, pars.speed, X=NULL, samples){
	smooth <- post$smooth
	estimpars <- post$estimpars
	
	N.samples <- dim(samples)[1]
	TT <- dim(samples)[3]
	xyrange <- c(-pars.speed$xmax/5, 1+pars.speed$xmax/5)
	T.past <- ncol(smooth$mu)

	ii.samples <- which(!is.na(samples[,1,1]))
	nn.sample <- length(ii.samples)
		
	cols.bg.sam <- viridis(TT, option='D', alpha=0.5)
	cols.bg.true <- viridis(TT, option='D')
	cols.fg.true <- rgb(0, 0.75, 1)
	
	# plot(estimpars$mu[1,], estimpars$mu[2,], xlim=xyrange, ylim=xyrange, pch='', xlab='x coordinate (m)', ylab='y coordinate (m)')
	plot(estimpars$mu[1,], estimpars$mu[2,], xlim=xyrange, ylim=xyrange, pch='', xlab='', ylab='', axes=F)
	scalebar2(0.2, 0, '20 cm', '', 'topleft')
	
	kk <- 1
	for (i in ii.samples){
		lines(samples[i,1,], samples[i,2,], t='o', cex=1, pch=16, col=viridis(nn.sample+1)[kk])
		kk <- kk + 1
	}

	if (!is.null(X)) {
		# lines(X[1,], X[2,], t='o', cex=1.3, pch=24, col=cols.fg.true, bg=cols.bg.true, lwd=2)
		lines(X[1,], X[2,], t='l', lwd=2, col=cols.fg.true)
		points(X[1,ii.samples], X[2,ii.samples], lwd=2, pch=24, col=cols.fg.true, bg=viridis(length(ii.samples)+1))		
	}
	
	# legend('topright', legend=c('samples', 'true trajectory'), pt.cex=c(1, 1.8), col=c(cols.bg.sam[round(TT/2)], cols.fg.true), pch=c(21, 24), pt.bg=c( cols.bg.sam[round(TT/2)], viridis(3)[2]), lwd=c(1,2))
}
