##########################################################
### speed functions
### 
##########################################################

## 1. smoothing
require(circular)
require(fmsb)

posfilter <- function(pos, sdfilt, graphics=F){
	## function to smooth the position data
	## INPUT: pos: matrix with the following columns: time (s); x pos (cm); ypos (cm), head direction (deg)
	## 	sdfilt: temporal filtering width (s)
	## 	graphics: T or F; whether reults should be plotted
	##
	## OUTPUT: txyvd, matrix with the following columns:
	## time (s); smoothed x pos (cm); smoothed ypos (cm), smoothed speed (cm/s); smoothed running direction (deg)
	
	# sdfilt <- 0.25 # s, the width of the Gaussian filter for position data
	dt_pos <- median(diff(pos[,1]))

	max_filt <- round(4 * sdfilt / dt_pos) * dt_pos
	xfilt <- seq(-max_filt, max_filt, by=dt_pos) # should have an odd length
	filt <- dnorm(xfilt, 0, sdfilt)
	filt <- filt / sum(filt)	
	if (graphics == T) plot(xfilt, filt, t='o')

	taxis <- seq(min(pos[,1], na.rm=T), max(pos[,1], na.rm=T), dt_pos)

	t.NA <- NULL
	ii_breaks <- which(diff(pos[,1]) > 10 * dt_pos)
	if (length(ii_breaks) > 0){
		tt_breaks <- matrix(NA, length(ii_breaks), 2)
		for (ibr in 1:length(ii_breaks)){
			tstart <- pos[ii_breaks, 1]
			tend <- pos[ii_breaks+1, 1]
			tt_breaks[ibr,] <- c(tstart, tend)
			istart <- min(which(taxis > tstart)) - floor(length(filt)/2)
			iend <- max(which(taxis < tend)) + floor(length(filt)/2)
			t.NA <- c(t.NA, istart:iend)
		}
	} else t.NA <- 0
	
	xxpos <- approx(pos[,1], pos[,2], taxis)$y
	xxpos <- c(rep(head(xxpos,1), floor(length(filt)/2)), xxpos, rep(tail(xxpos,1), floor(length(filt)/2)))
	xpos <- filter(xxpos, filt)
	xpos <- xpos[!is.na(xpos)]
	xpos[t.NA] <- NA
	
	yypos <- approx(pos[,1], pos[,3], taxis)$y
	yypos <- c(rep(head(yypos,1), floor(length(filt)/2)), yypos, rep(tail(yypos,1), floor(length(filt)/2)))
	ypos <- filter(yypos, filt)
	ypos <- ypos[!is.na(ypos)]
	ypos[t.NA] <- NA
	
	if (graphics == T){
		speed <- sqrt(diff(pos[,2])^2 + diff(pos[,3])^2) / diff(pos[,1])
		speed <- c(speed[1], speed)
	}
	
	motion_dir <- atan2(diff(ypos), diff(xpos))
	motion_dir <- c(motion_dir[1], motion_dir)
	
	speed2 <- sqrt(diff(xpos)^2 + diff(ypos)^2) / dt_pos
	speed2 <- c(speed2[1], speed2)
	
	txyvd <- cbind(taxis, xpos , ypos, speed2, motion_dir)
	
	if (graphics == T){
		tstart <- floor(min(taxis))
		tend <- tstart + 300
		jstart <- min(which(taxis > tstart))
		jend <- max(which(taxis < tend))
		jj <- jstart:jend

		istart <- min(which(pos[,1] > tstart))
		iend <- max(which(pos[,1] < tend))
		ii <- istart:iend

			
		par(mfcol=c(3,1)); par(mar=c(2,4,2,2))
		plot(pos[ii,1], pos[ii,2], t='l', axes=F, xlab='', ylab='', ylim=c(0, 200)); axis(2, las=2)
		lines(taxis[jj], xpos[jj], col=3)
		
		lines(pos[ii,1], pos[ii,3], t='l', lty=3)
		lines(taxis[jj], ypos[jj], col=4, lty=1)
		
		par(mar=c(4,4,2,2))
		plot(pos[ii,1], speed[ii], t='l', axes=F, xlab='', ylab='', ylim=c(0, max(speed, na.rm=T))); axis(2, las=2); axis(1)
		lines(taxis[jj], speed2[jj], col=3)
		
		plot(xpos[jj], ypos[jj], t='l', col=3)
		lines(pos[ii,2], pos[ii,3], t='l', col=1)
	}

	txyvd
}

# T_theta_min <- 1
# T_SPW_min <- 0.25

get_long_periods <- function(states, N_true_min=29.5, N_false_min=7.5){
	# a function to detect periods of TRUEs and FALSEs
	# minimum duration of TRUEs: N_true_min
	# minimum interval between TRUE states: N_false_min

	# first, remove short FALSE periods by making them TRUE
	
	states[is.na(states)] <- F
	
	nStates <- rle(states)
	iStates <- cumsum(nStates$lengths)

	for (j in 1:length(iStates)){
		if ((nStates$value[j] == F) & (nStates$lengths[j] <= N_false_min)) {# FALSE, and short
			iend <- iStates[j]
			istart <- iend - nStates$lengths[j] + 1
			states[istart:iend] <- TRUE
		}
	}

	# second, remove short TRUE periods by making them FALSE
	nStates <- rle(states)
	iStates <- cumsum(nStates$lengths)

	for (j in 1:length(iStates)){
		if ((nStates$value[j] == T) & (nStates$lengths[j] <= N_true_min)) {# not hteta, but short
			iend <- iStates[j]
			istart <- iend - nStates$lengths[j] + 1
			states[istart:iend] <- FALSE
		}
	}
	
	states
	
}

##########################################################
## naive likelihood - iid samples from the prior 
## - to make it comparable to the static and dynamic, we need to calculate it as if it was a continuous probability density P(x,y)
## the support of the PD is the 200x200 cm arena, so when we take the integral \int_0^200 \int_0^200 P(x,y) dx dy = 1
## we approximate it with a step-like distribution, where the P(x,y) is constant within 5 cm x 5 cm bins
## this is the reason for the division by 25 in P_map <- trainmap / sum(trainmap) / 25

LL_pos <- function(xy=NULL, testmap=NULL, trainmap){
	LL <- 0
	P_map <- trainmap / sum(trainmap) / 25	
	
	dt_pos <- attr(trainmap, 'dt')
	if (is.null(testmap)){
		if (is.null(xy)) stop('testmap or xy should be given!')
		xw <- diff(attr(trainmap, 'xcenters')[1:2])
		yw <- diff(attr(trainmap, 'ycenters')[1:2])
		ix_run <- xy[,1] %/% xw + 1
		iy_run <- xy[,2] %/% yw + 1
		for (i in 1:length(ix_run)){
			LL <- LL + log(P_map[ix_run[i], iy_run[i]])
		}	
	} else {
		for (i in 1:nrow(testmap)){
			for (j in 1:ncol(testmap)){
				LL <- LL + log(P_map[i, j]) * testmap[i,j] / dt_pos
			}
		}
	}
	LL
}


##########################################################
## plotting motion statistics ...
ExpFitErr <- function(tau, x, ref){
	expfun <- exp(-x/tau)
	err <- sum((expfun - ref)^2)
	err
}

plot_motion_stats <- function(xyvd, dt_pos, file=NULL){
	if (!is.null(file)){
		pdf(file=file, 10, 5, useD=F)
	}
	par(mfrow=c(2,4))
	col1 <- viridis(5, option='B')[4]
	col2 <- viridis(5, option='D')[4]

	vvr <- xyvd[,3] # velocity during running
	mvvr <- mean(vvr, na.rm=T) # mean velocity during running
	
	max_brs <- max(100, ceiling(max(vvr, na.rm=T)/2)*2)
	brs <-seq(0, max_brs, by=2)
	hist(vvr, br=brs, freq=F, main=paste('speed, mean: ', round(mvvr), 'cm/s'), xlab='speed (cm/s)')
	abline(v=mvvr, col=col1, lwd=3)
	
	acf_vvr <- acf(vvr, na.action=na.pass, lag=200, plot=F)
	tau <- optim(par=1, ExpFitErr, gr=NULL, x=acf_vvr$lag*dt_pos, ref=acf_vvr$acf, method='BFGS')$par
	plot(acf_vvr$lag*dt_pos, acf_vvr$acf, t='l', lty=1, lwd=2, axes=F, xlab='time (s)', ylab='speed autocorrelation', main=paste('Tau =', round(tau, 3), 's'), ylim=c(0, 1))
	lines(acf_vvr$lag*dt_pos, exp(-acf_vvr$lag*dt_pos/tau), col=col1, lwd=2)
	axis(1); axis(2, las=2)
	
	###############################################
	### speed changes - only during running
	
	dvvr <- diff(vvr) / dt_pos
	mdvvr <- mean(dvvr, na.rm=T) # mean acceleration - 0
	vdvvr <- mean(dvvr^2, na.rm=T) - mean(dvvr, na.rm=T)^2 # SD of accelearation
	
	
	max_brs <- max(120, ceiling(max(abs(dvvr), na.rm=T)/10)*10)
	brs <-seq(-max_brs, max_brs, by=5)
	hist(dvvr, br=brs, freq=F, main=paste('acceleration; SD=', round(sqrt(vdvvr)), 'cm/s^2'), xlab='acceleration (cm / s^2)')
	lines(brs, dnorm(brs, mdvvr, sqrt(vdvvr)), col=col1, lwd=2)
	
	acf_a <- acf(dvvr, na.action=na.pass, lag=200, plot=F)
	tau_a <- optim(par=1, ExpFitErr, gr=NULL, x=acf_a$lag*dt_pos, ref=acf_a$acf, method='BFGS')$par
	plot(acf_a$lag*dt_pos, acf_a$acf, t='l', lty=1, lwd=2, axes=F, xlab='time (s)', ylab='HD autocorrelation', main=paste('Tau =', round(tau_a, 3), 's'))
	lines(acf_a$lag*dt_pos, exp(-acf_a$lag*dt_pos/tau_a), col=col1, lwd=2)
	axis(1); axis(2, las=2)
	
	###############################################
	## head direction, and its acf
	hd <- xyvd[,4]
	brs <- seq(-pi-pi/16, pi+pi/16, by=pi/8)
	h_hd <- hist(hd, br=brs, plot=F)
	data <- as.data.frame(matrix(h_hd$density[1:16], ncol=16))
	data[1] <- data[1] + h_hd$density[17]
	colnames(data) <- c('-pi', '', '', '', '-pi/2', '', '', '', '0', '', '', '', 'pi/2', '', '', '')
	data <- rbind(rep(0.25, 16), rep(0, 16), data)
	radarchart(data, axistype=1, pcol=viridis(5)[3], pfcol=viridis(5, alpha=0.5)[4], plwd=4, cglcol="grey", cglty=1, axislabcol="grey", caxislabels=0:4/4, cglwd=0.8,  vlcex=0.8)
	mtext('head direction', 3, 1.5, font=2, cex=0.8)
	
	acf_hd <- acf(hd, na.action=na.pass, lag=200, plot=F)
	tau_hd <- optim(par=1, ExpFitErr, gr=NULL, x=acf_hd$lag*dt_pos, ref=acf_hd$acf, method='BFGS')$par
	plot(acf_hd$lag*dt_pos, acf_hd$acf, t='l', lty=1, lwd=2, axes=F, xlab='time (s)', ylab='HD autocorrelation', main=paste('Tau =', round(tau_hd, 3), 's'))
	lines(acf_hd$lag*dt_pos, exp(-acf_hd$lag*dt_pos/tau_hd), col=col2, lwd=2)
	axis(1); axis(2, las=2)
	
	###############################################
	### direction changes - only during running
	dir_change <- diff(hd)
	ii_small <- which(dir_change < (-1) * pi)
	ii_big <- which(dir_change > pi)
	dir_change[ii_small] <- dir_change[ii_small] + 2*pi
	dir_change[ii_big] <- dir_change[ii_big] - 2*pi
	
	dcr <- dir_change
	
	brs <- seq(-3.2, 3.2, by=0.02)
	md <- mean(dcr, na.rm=T)
	vd <- mean(dcr^2, na.rm=T) - mean(dcr, na.rm=T)^2
	hist(dcr, br=brs, freq=F, xlim=c(-3*sqrt(vd), 3*sqrt(vd)), main='heading change', xlab='heading change (rads)')
	lines(brs, dnorm(brs, md, sqrt(vd)), col=col2, lwd=2)
	
	
	acf_hdc <- acf(dcr, na.action=na.pass, lag=200, plot=F)
	tau_hdc <- optim(par=1, ExpFitErr, gr=NULL, x=acf_hdc$lag*dt_pos, ref=acf_hdc$acf, method='BFGS')$par
	plot(acf_hdc$lag*dt_pos, acf_hdc$acf, t='l', lty=1, lwd=2, axes=F, xlab='time (s)', ylab='HD change autocorrelation', main=paste('Tau =', round(tau_hdc, 3), 's'))
	lines(acf_hdc$lag*dt_pos, exp(-acf_hdc$lag*dt_pos/tau_hdc), col=col2, lwd=2)
	axis(1); axis(2, las=2)
	
	if (!is.null(file)){
		dev.off()
	}
	
	stats <- list(mu_v= mvvr, mu_a=mdvvr, sd_a=sqrt(vdvvr))
}


##########################################################
## static likelihood - previous position + Gaussian random walk
## xy: position matrix with columns x and y
## i_data: index of datapoints used for evaluation
## sd_xy: SD of the Gaussian random walk - mean distance between neighbouring points
## trainmap: occupancy map of the animal during training time

LL_static <- function(xy, i_data, sd_xy, trainmap){
	# first point from the prior
	P_map <- trainmap / sum(trainmap) / 25
	i_start_end <- find_chunks(i_data)
	N.runs <- nrow(i_start_end)
	dt_pos <- attr(trainmap, 'dt')

	xw <- diff(attr(trainmap, 'xcenters')[1:2])
	yw <- diff(attr(trainmap, 'ycenters')[1:2])
	const <- (-1) * log(2*pi*(sd_xy* dt_pos)^2) # 1/sqrt( (2*pi*sd_xy^2)^2)
	icov <- diag(1/(sd_xy*dt_pos)^2, 2)

	LL <- 0
	
	for (j in 1:N.runs){
		xy.j <- xy[i_start_end[j,1]:i_start_end[j,2],]
		# first point
		xy.past <- xy.j[1,]
		ix_run <- xy.past[1] %/% xw + 1
		iy_run <- xy.past[2] %/% yw + 1
		LL <- LL + log(P_map[ix_run, iy_run])
	
		# all other points
		N.j <- nrow(xy.j)
		LL <- LL + (N.j-1) * const
		for (i in 2:N.j){
			xyt <- xy.j[i,]
			LL <- LL - 1/2 * t(xyt - xy.past) %*% icov %*% (xyt - xy.past)
			xy.past <- xyt
		}		
	}
	LL
}

##########################################################
## dynamic likelihood - previous position and SPEED + Gaussian random walk on SPEED
## we assume that the motion direction is kept constant
## xy: position matrix with columns x, y
## i_data: index of datapoints used for evaluation
## sd_xy: SD of the Gaussian random walk - mean distance between neighbouring points - used to predict the second point in the trajectory
## sd_a: SD of the Gaussian random walk of the SPEED (we assume the same noise on the parallel and the orthogonal speed components)
## trainmap: occupancy map of the animal during training time

LL_dynamic <- function(xy, i_data, sd_xy, sd_a, trainmap){
	# first point from the prior
	P_map <- trainmap / sum(trainmap) / 25
	i_start_end <- find_chunks(i_data)
	N.runs <- nrow(i_start_end)
	dt_pos <- attr(trainmap, 'dt')

	xw <- diff(attr(trainmap, 'xcenters')[1:2])
	yw <- diff(attr(trainmap, 'ycenters')[1:2])
	LL <- 0

	const_xy <- (-1) * log(2*pi*(sd_xy*dt_pos)^2) # 1/sqrt( (2*pi*sd_xy^2)^2)
	icov_xy <- diag(1/(sd_xy*dt_pos)^2, 2)

	sd_pos <- sd_a * dt_pos^2 # next = past + v*dt; v = v0 + dv/dt * dt = v0 + a*dt; SD(next) = SD(a) * dt^2
	const_pos <- (-1) * log(2*pi*sd_pos^2) # 1/sqrt( (2*pi*sd_a^2)^2)
	icov_pos <- diag(1/sd_pos^2, 2)

	for (j in 1:N.runs){
		xy.j <- xy[i_start_end[j,1]:i_start_end[j,2],]

		## first point - iid
		xy.past <- xy.j[1,]
		ix_run <- xy.past[1] %/% xw + 1
		iy_run <- xy.past[2] %/% yw + 1
		LL <- c(LL, log(P_map[ix_run, iy_run]))
	
		# second point - from static
		xy.recent <- xy.j[2,]
		LL <- c(LL, - 1/2 * t(xy.recent - xy.past) %*% icov_xy %*% (xy.recent - xy.past) + const_xy)
		
		# all other points point - from dynamic
		N.j <- nrow(xy.j)
		# LL <- LL + (N.j-2) * const_pos

		for (i in 3:N.j){
			# xy = xy + v * dt where the velocity is random, approximately Gaussian
			
			xyt <- xy.j[i,]
			v_past <- xy.recent - xy.past
			mu_next <- xy.recent + v_past
			LL <- c(LL,  - 1/2 * t(xyt - mu_next) %*% icov_pos %*% (xyt - mu_next) + const_pos)
			xy.past <- xy.recent
			xy.recent <- xyt
		}
	}
	LL
}

