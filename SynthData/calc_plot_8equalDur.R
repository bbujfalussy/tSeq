## calculate the (forward and lateral) bias of the decoder and the firing rate as a function of theta phase
## - using standard ratemaps
## - 8 temporal windows of 120deg width each

require(ellipse)
ell2 <- ellipse::ellipse

dir.create('./SimFigs/Fig3_4/', showW=F)

predcos <- function(pars, xx){
	amp <- pars[1]
	baseline <- pars[2]
	phase <- pars[3]
	pred <- baseline + amp * cos(xx-phase)
	pred	
}

errcos <- function(pars, xx, yy){
	pred <- predcos(pars, xx)
	err <- sum((pred - yy)^2)
	err
}

cat('calculating theta chunks...')
theta_phases <- seq(0, 7/8, by=1/8)
phase_ends <- theta_phases + 1/3

theta_starts <- theta.starts.orig[1:N.theta]
theta_ends <- theta.starts.orig[-1]
bin_lengths <- theta_ends - theta_starts


bin.start.end1 <- cbind(theta_starts, theta_starts + 1/3* bin_lengths)
rast1 <- spt2raster_bins(spt, bin.start.end1, N.cells=N.cells)

bin.start.end2 <- cbind(theta_starts + 1/8 * bin_lengths, theta_starts + (1/3+1/8)* bin_lengths)
rast2 <- spt2raster_bins(spt, bin.start.end2, N.cells=N.cells)

bin.start.end3 <- cbind(theta_starts + 2/8 * bin_lengths, theta_starts + (1/3+2/8)* bin_lengths)
rast3 <- spt2raster_bins(spt, bin.start.end3, N.cells=N.cells)

bin.start.end4 <- cbind(theta_starts + 3/8 * bin_lengths, theta_starts + (1/3+3/8)* bin_lengths)
rast4 <- spt2raster_bins(spt, bin.start.end4, N.cells=N.cells)

bin.start.end5 <- cbind(theta_starts + 4/8 * bin_lengths, theta_starts + (1/3+4/8)* bin_lengths)
rast5 <- spt2raster_bins(spt, bin.start.end5, N.cells=N.cells)

bin.start.end6 <- cbind(theta_starts + 5/8 * bin_lengths, theta_starts + (1/3+5/8)* bin_lengths)
rast6 <- spt2raster_bins(spt, bin.start.end6, N.cells=N.cells)

bin.start.end7 <- cbind(theta_starts[1:N.theta-1] + 6/8 * bin_lengths[1:N.theta-1], theta_starts[2:N.theta] + (2/24)* bin_lengths[2:N.theta])
rast7 <- spt2raster_bins(spt, bin.start.end7, N.cells=N.cells)

bin.start.end8 <- cbind(theta_starts[1:N.theta-1] + 7/8 * bin_lengths[1:N.theta-1], theta_starts[2:N.theta] + (5/24)* bin_lengths[2:N.theta])
rast8 <- spt2raster_bins(spt, bin.start.end8, N.cells=N.cells)

estimfile <- paste('estimpos8_', code, '_spPth', spPth, reg, '_', jitter, '.RData', sep='')
if (estimfile %in% list.files('./SimFigs/Fig3_4/')){
	estimf <- paste('./SimFigs/Fig3_4/', estimfile, sep='')	
	load(file=estimf)
	cat('loading decoded position ...')
	estimpos1 <- estimp$est1
	estimpos2 <- estimp$est2
	estimpos3 <- estimp$est3
	estimpos4 <- estimp$est4
	estimpos5 <- estimp$est5
	estimpos6 <- estimp$est6
	estimpos7 <- estimp$est7
	estimpos8 <- estimp$est8
	
	estimp <- list(est1=estimpos1, est2=estimpos2, est3=estimpos3, est4=estimpos4, est5=estimpos5, est6=estimpos6, est7=estimpos7, est8= estimpos8)
} else {
	cat('decoding position ...')
	estimpos1 <- decode_ratemaps(rast1, bin.start.end1[,2] - bin.start.end1[,1], ratemaps, posmap)
	estimpos2 <- decode_ratemaps(rast2, bin.start.end2[,2] - bin.start.end2[,1], ratemaps, posmap)
	estimpos3 <- decode_ratemaps(rast3, bin.start.end3[,2] - bin.start.end3[,1], ratemaps, posmap)
	estimpos4 <- decode_ratemaps(rast4, bin.start.end4[,2] - bin.start.end4[,1], ratemaps, posmap)
	estimpos5 <- decode_ratemaps(rast5, bin.start.end5[,2] - bin.start.end5[,1], ratemaps, posmap)
	estimpos6 <- decode_ratemaps(rast6, bin.start.end6[,2] - bin.start.end6[,1], ratemaps, posmap)
	estimpos7 <- decode_ratemaps(rast7, bin.start.end7[,2] - bin.start.end7[,1], ratemaps, posmap)
	estimpos8 <- decode_ratemaps(rast8, bin.start.end8[,2] - bin.start.end8[,1], ratemaps, posmap)
	
	estimp <- list(est1=estimpos1, est2=estimpos2, est3=estimpos3, est4=estimpos4, est5=estimpos5, est6=estimpos6, est7=estimpos7, est8= estimpos8)
	
	estimf <- paste('./SimFigs/Fig3_4/', estimfile, sep='')	
	save(estimp, file=estimf)	
}

motion_dir <- atan2(diff(txy[,3]), diff(txy[,2]))
L <- nrow(txy)
t_motion_dir <- (txy[1:(L-1),1] + txy[2:L]) / 2

if (graphics){
	cols <- colormap(colormap=colormaps$hsv, 8, alpha=0.1)
	mcols <- colormap(colormap=colormaps$hsv, 18)[c(1, 3, 4, 5, 9, 11, 13, 15)]

	##################################################
	## aligned decoded positions, number of spikes and bin duration
	qq <- 0 #19/20 #19/20
	aligned_errs <- list()
	max_nsp <- 0
	
	for (j in 1:8){
		
		if (j == 1) {estimpos <- estimpos1; rast <- rast1; bin.start.end <- bin.start.end1}
		if (j == 2) {estimpos <- estimpos2; rast <- rast2; bin.start.end <- bin.start.end2}
		if (j == 3) {estimpos <- estimpos3; rast <- rast3; bin.start.end <- bin.start.end3}
		if (j == 4) {estimpos <- estimpos4; rast <- rast4; bin.start.end <- bin.start.end4}
		if (j == 5) {estimpos <- estimpos5; rast <- rast5; bin.start.end <- bin.start.end5}
		if (j == 6) {estimpos <- estimpos6; rast <- rast6; bin.start.end <- bin.start.end6}
		if (j == 7) {estimpos <- estimpos7; rast <- rast7; bin.start.end <- bin.start.end7}
		if (j == 8) {estimpos <- estimpos8; rast <- rast8; bin.start.end <- bin.start.end8}
		
		nsp <- rowSums(rast)
		# rates <- nsp / (bin.start.end[,2]-bin.start.end[,1]) / N.cells
		t_bin <- bin.start.end[,2]-bin.start.end[,1]
		
		refx <- approx(txy[,1], txy[,2], rowMeans(bin.start.end))$y / 100
		refy <- approx(txy[,1], txy[,3], rowMeans(bin.start.end))$y / 100
		ref_dir <- approx(t_motion_dir, motion_dir, rowMeans(bin.start.end), rule=2)$y
		
		abs_estim_dir <- atan2(estimpos$mean[,2] - refy, estimpos$mean[,1] - refx)
		estim_dir <- (abs_estim_dir - ref_dir + pi) %% (2*pi) - pi										
		rsqerror <- sqrt(rowSums((estimpos$mean - cbind(refx, refy))^2))

		xi <- rsqerror * cos(estim_dir)
		yi <- rsqerror * sin(estim_dir)
		xynt <- cbind(xi, yi, nsp, t_bin)
		if (max(xynt[,3]) > max_nsp) max_nsp <- max(xynt[,3])
		aligned_errs[[j]] <- xynt
	}

	##################################################
	## plot the decoded positions and thear mean and spread
	brs <- 0 :(max_nsp+1) - 0.5
	hnsp <- matrix(NA, 8, max_nsp+1)
	for (j in 1:8){
		hnsp[j,] <- hist(aligned_errs[[j]][,3], br=brs, plot=F)$counts
	}
	image(1:8, 0:max_nsp, hnsp, axes=F, ylab='', xlab='')
	lines(1:8, hnsp %*% 0: max_nsp / rowSums(hnsp), lwd=2)
	axis(2, las=2)
	axis(1, 1:8, 1:8*45)
	min_counts <- apply(hnsp, 2, min)
	plot(0:max_nsp, cumsum(min_counts) / sum(min_counts))
	spike_threshold <- min(which(cumsum(min_counts) / sum(min_counts) > qq)) - 1
	abline(h=qq, v=spike_threshold)
	
	## function to perform thinning
	thin_vec <- function(nsp, counts_thin){
		# nsp: spike counts in each bin to downsample
		# counts_thin: histogram, the number of elements to keep from 0 to Nmax
		Nmax <- length(counts_thin) - 1
		M <- length(counts_thin)
		index_keep <- NA
		for (j in 1:M){
			if (counts_thin[j] > 0){
				index_j <- which(nsp==(j-1))
				index_keep <- c(index_keep, sample(index_j, counts_thin[j]))
			}
		}
		sort(index_keep)
	}
	
	shifter <- function(x, n = 1) {
	  if (n == 0) x else c(tail(x, -n), head(x, n))
	}
	
	#############################################################################################
	## to make a fair comparison, we need to perform thinning - downsample nspike-histogram to equaliye counts accross phase
	
	ms <- matrix(NA, 8, 2)
	spread <- rep(NA, 8)
	rate <- matrix(NA, 8, 2)
	
	spike_threshold <- min(which(cumsum(min_counts) / sum(min_counts) > qq)) - 1
	
	filename <- paste('./SimFigs/Fig3_4/postmean_phase_', code, '_spPth', spPth, reg, '_', jitter, '_q', qq, '.pdf', sep='')
	pdf(file=filename, 4, 4, useD=F)
	par(mar=c(1,1,1,1))
	xylist <- list()
	plot(0, 0, xlim=c(-0.3, 0.3), ylim=c(-0.3,0.3), t='n', xlab='', ylab='', axes=F)
	
	for (j in 1:8){
		xynt <- aligned_errs[[j]]
		all_nsp <- xynt[,3]
		all_t_bin <- xynt[,4]
		jj <- thin_vec(all_nsp, min_counts)
		xy <- xynt[jj,1:2]
		nsp <- xynt[jj,3]
		ii <- which(nsp > spike_threshold)
		xylist[[j]] <- xy[ii,]
		ms[j,] <- colMeans(xy[ii,])
		spread[j] <- det(cov(xy[ii,]))^(1/4)
		rate[j,1] <- sum(all_nsp) / sum(all_t_bin) / N.cells	
		rate[j,2] <- sd(all_nsp / all_t_bin / N.cells)
		points(xy[ii,1], xy[ii,2], pch=16, cex=1, col=cols[j])
	}
	for (j in 1:8){
		lines(ell2(cov(xylist[[j]]), centre=c(ms[j,1], ms[j,2]), level=0.5), t='l', col=mcols[j], lwd=2)
		points(ms[j,1], ms[j,2], pch=21, col=1, bg=mcols[j], cex=2)	
	}

	scalebar2(0.05, 0.05, ' ', '5 cm')
	points(0, 0, pch=3, cex=3, lwd=2)
	dev.off()
		

#########################################
## bias, spread and number of spikes as a function of theta phase
	ms_cm <- ms * 100	
	spread_cm <- spread * 100	

	initpars <- c(6, 6, 0)
	theta_phases	<- seq(0, by=2*pi/8, length=8)
	newpars <- optim(initpars, errcos, gr=NULL, xx=theta_phases, yy=ms_cm[,1])
	xx <- seq(0, 2*pi, length=100)
	yy <- predcos(newpars$par, xx)


	pchs <- c(1, 0, 5, 2)	
	peakPhase <- newpars$par[3] %% (2*pi)
	if (newpars$par[1] < 0) peakPhase <- round((peakPhase + pi) %% (2*pi), 2) else peakPhase <- round(peakPhase, 2)

	filename <- paste('./SimFigs/Fig3_4/bias_rateM_', code, '_spPth', spPth, reg, '_', jitter, '_q', qq, '.pdf', sep='')
	pdf(file=filename, useD=F, 2.5, 5)
	par(mfcol=c(3,1))
	par(mar=c(1, 4, 4, 1))
	
	captiontext <- paste('peak of phase:', peakPhase, 'rad')
	plot(c(theta_phases, 2*pi), c(ms_cm[,1], ms_cm[1,1]), col=1, pch=20+i.code, t='o', bg=mcols, xlab='', ylab='decoding bias (cm)', axes=F, ylim=c(-20, 10), xlim=c(0, 2*pi), main=captiontext, cex=1.5)
	# points(c(theta_phases, 2*pi), c(ms_cm[,2], ms_cm[1,2]), col=1, pch=22, bg=mcols, t='o')
	# lines(xx, yy)
	axis(2, las=2)
	legend('topright', leg=c('forward', 'lateral'), pch=c(21, 22), pt.bg=mcols[2], bty='n')
	
	plot(c(theta_phases, 2*pi), c(spread_cm, spread_cm[1]), col=1, pch=20 + i.code, t='o', bg=mcols, xlab='', ylab='decoding spread (cm)', axes=F, xlim=c(0, 2*pi), ylim=c(8, 18), main='120 deg window', cex=1.5)
	axis(2, las=2)
	
	
	par(mar=c(4, 4, 1, 1))
	# plotCI(c(theta_phases, 2*pi), c(rate[,1], rate[1,1]), c(rate[,2], rate[1,2]), gap=0, ylim=c(0, 3), axes=F, xlab='theta phase (deg)', ylab='firing rate (Hz)', pch=pchs[i.code])
	plot(c(theta_phases, 2*pi), c(rate[,1], rate[1,1]), ylim=c(0, 3), axes=F, xlab='theta phase (deg)', ylab='firing rate (Hz)', pch=20 + i.code, t='o', bg=mcols, cex=1.5)
	axis(1, c(0, pi, 2*pi), c(0, 180, 360))
	axis(2, 0:3, las=2)
	# points(c(theta_phases, 2*pi), c(rate[,1], rate[1,1]), pch=20 + i.code, bg=mcols, t='o')
	yy2 <- 3 * (yy - min(yy)) / diff(range(yy)) 
	lines(xx, yy2)
	dev.off()
	
	simdata <- list(ms=ms_cm, spread=spread_cm, rate=rate)
	filename <- paste('./SimFigs/Fig3_4/bias_rateM_', code, '_spPth', spPth, reg, '_', jitter, '_q', qq, '.RData', sep='')
	save(simdata, file=filename)
	
}