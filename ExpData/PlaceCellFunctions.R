### place cell functions
library(viridis)
library(spatstat)
library(gtools)
library(colormap)
library("Matrix")
source('~/Programs/Code/HPC_Data/SimRunND.R', chdir = TRUE)

# spi <- spikes[spikes[,2] == 1,1]
# txyv <- cbind(pos[,1], xpos , ypos, speed2)
# vmin <- 5
# xw <- 10; yw <- 10
# sigma <- 10

estimate_ratemap <- function(spi, txy, i_data, xbrs, ybrs, sigma, posmap, graphics=F, dim=2){
## a function to estimate the ratemap from spikes during exploration
## spikes: matrix with the spike times (sec) and cell ids as columns
## pos: matrix with time (sec), smoothed x and y position (cm) as columns 
## i_data: index for the position to include in the estimate
## sigma: width of a Gaussian smoothing kernel
## xw, yw: bin width for counting spikes 2D histogram; arena is 200 x 200 cm

	ii_run <- approx(txy[,1], i_data, spi, method='constant')$y
	spi_run <- spi[!!ii_run]

	xw <- median(diff(xbrs))
	yw <- median(diff(ybrs))
	
	if (xw - max(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
	if (xw - min(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
	if (yw - max(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')
	if (yw - min(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')

	txy2 <- txy
	txy2[,2] <- txy[,2] - xbrs[1]
	x_spikes <- approx(txy2[,1], txy2[,2], spi_run)$y
	ix_spikes <- x_spikes %/% xw + 1
	if (dim == 1){
		Nbin <- length(xbrs) - 1
		countmap <- matrix(0, 1, Nbin)
		for (i in 1: length(ix_spikes)){
			countmap[ix_spikes[i]] <- countmap[ix_spikes[i]] + 1
		}

		cm2 <- rbind(countmap, countmap)
		rm2 <- blur(as.im(cm2), sigma/xw, normalise=T)$v
		ratemap <- rm2[1,] / posmap
		if (graphics){
			plot(xbrs[-1] - xw/2, countmap/posmap, xlab='x coordinate', ylab='density', t='o', pch=15)
			lines(xbrs[-1] - xw/2, ratemap, col=viridis(5)[4], lwd=2)
			points(spi_run, runif(length(spi_run), 0, max(ratemap)), pch=16, cex=0.5, col=grey(0.8, alpha=0.5))
		}

	} else {
		txy2[,3] <- txy[,3] - ybrs[1]
		y_spikes <- approx(txy2[,1], txy2[,3], spi_run)$y
		iy_spikes <- y_spikes %/% yw + 1
		
		Nbin <- c(length(xbrs), length(ybrs)) - 1
		countmap <- matrix(0, Nbin, Nbin)
		for (i in 1: length(ix_spikes)){
			countmap[ix_spikes[i], iy_spikes[i]] <- countmap[ix_spikes[i], iy_spikes[i]] + 1
		}
		
		# if (is.null(posmap)) posmap <- estimate_posmap(txy, vv, vmin, xw, yw, sigma)
		ratemap <- blur(as.im(countmap), sigma/xw, normalise=T)$v / posmap
	
		if (graphics>0){
			par(mar=c(1,1,1,1)/4)
			image(xbrs[-1] - xw/2, ybrs[-1]-yw/2, ratemap, col=colormap(colormaps$hot, 48), xlab='', ylab='', axes=F)
			# image(1:20*10-5, 1:20*10-5, ratemap, col=colormap(colormaps$hot, 48), xlab='', ylab='', axes=F)
			# image(1:20*10-5, 1:20*10-5, poscount, col=viridis(24))
			if (graphics>1) points(x_spikes + xbrs[1], y_spikes + ybrs[1], pch=16, cex=0.4, col=viridis(6, option='D', alpha=0.5)[4])
		}
	}
	
	attr(ratemap, 'nsp') <- length(x_spikes)
	attr(ratemap, 'xcenters') <- xbrs[-1] - xw/2
	if (dim==2) attr(ratemap, 'ycenters') <- ybrs[-1] - yw/2
	ratemap
}

estimate_ratemap_xyw <- function(counts, xyw, xbrs, ybrs, sigma, posmap, graphics=F, dim=2){
## a function to estimate the ratemap from spikes during exploration from theta phases
## spikes: matrix with the spike times (sec) and cell ids as columns
## pos: matrix with time (sec), smoothed x and y position (cm) as columns 
## i_data: index for the position to include in the estimate
## sigma: width of a Gaussian smoothing kernel
## xw, yw: bin width for counting spikes 2D histogram; arena is 200 x 200 cm
	
	if (dim==1){
		xw <- median(diff(xbrs))
		
		if (xw - max(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
		if (xw - min(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
	
		xyw2 <- xyw
		xyw2[,1] <- xyw[,1] - xbrs[1]
	
		i_spikes <- which(counts > 0)
		spikes <- counts[i_spikes]
		xyw_sp <- xyw2[i_spikes,]
		
		ix <- xyw_sp[,1] %/% xw + 1
	
		Nbin <- length(xbrs) - 1
		countmap <- rep(0, Nbin)
		for (i in 1: length(ix)){
			countmap[ix[i]] <- countmap[ix[i]] + spikes[i]
		}
		
		cm2 <- rbind(countmap, countmap)
		rm2 <- blur(as.im(cm2), sigma/xw, normalise=T)$v
		ratemap <- rm2[1,] / posmap
	
		if (graphics>0){
			par(mar=c(1,1,1,1)/4)
			plot(xbrs[-1] - xw/2, ratemap, xlab='', ylab='', axes=F)
		}
		attr(ratemap, 'nsp') <- sum(counts)
		attr(ratemap, 'xcenters') <- xbrs[-1] - xw/2
	} else {
		xw <- median(diff(xbrs))
		yw <- median(diff(ybrs))
		
		if (xw - max(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
		if (xw - min(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
		if (yw - max(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')
		if (yw - min(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')
	
		xyw2 <- xyw
		xyw2[,1] <- xyw[,1] - xbrs[1]
		xyw2[,2] <- xyw[,2] - ybrs[1]
	
		i_spikes <- which(counts > 0)
		spikes <- counts[i_spikes]
		xyw_sp <- xyw2[i_spikes,]
		
		ix <- xyw_sp[,1] %/% xw + 1
		iy <- xyw_sp[,2] %/% yw + 1
	
		Nbin <- c(length(xbrs), length(ybrs)) - 1
		countmap <- matrix(0, Nbin[1], Nbin[2])
		for (i in 1: length(ix)){
			countmap[ix[i], iy[i]] <- countmap[ix[i], iy[i]] + spikes[i]
		}
		
		# if (is.null(posmap)) posmap <- estimate_posmap(txy, vv, vmin, xw, yw, sigma)
		ratemap <- blur(as.im(countmap), sigma/xw, normalise=T)$v / posmap
	
		if (graphics>0){
			par(mar=c(1,1,1,1)/4)
			image(xbrs[-1] - xw/2, ybrs[-1]-yw/2, ratemap, col=colormap(colormaps$hot, 48), xlab='', ylab='', axes=F)
		}
		
		attr(ratemap, 'nsp') <- sum(counts)
		attr(ratemap, 'xcenters') <- xbrs[-1] - xw/2
		attr(ratemap, 'ycenters') <- ybrs[-1] - yw/2
	}
	ratemap
}


## a function to plot a few statistics about the place cells
## intputs: ratemaps: array of dimension N x X x Y with the X x Y ratemap matrix of the N neurons
## Ecells index of the Excitatory cells - putative place cells
## errFisher: error estimate calculated from Fisher information, in m
## cellids: the id of cells to plot the place fields
## optional  - nspikes: number of spikes for the cells during the sampled interval
## Ttotal: the duration of the sampled interval

plot_rate_stats <- function(ratemaps, actEcells, Ecells, cellids, nspikes=NULL, Ttotal=NULL, posmap=NULL){
	par(mfrow=c(2,4))
	par(mar=c(4, 4, 2, 1))
	ncells <- dim(ratemaps)[1]
	
	PCCol <- viridis(6, alpha=1)[3] # green
	aPCCol <- viridis(6, alpha=1)[5] # green
	IntCol <- grey(0.7)
	iCol <- viridis(6, alpha=1, option='B')[5] # orange

	dx <- attr(ratemaps, 'xcenters')[2]-attr(ratemaps, 'xcenters')[1]
	if (is.null(nspikes)) {
		meanrate <- mean(apply(ratemaps[Ecells,,], 1, mean))
		rate1 <- apply(ratemaps, 1, mean)
	} else {
		meanrate <- sum(nspikes[Ecells]) / Ttotal / length(Ecells)
		rate1 <- nspikes / Ttotal
	}
	title <- paste('r=', round(meanrate, 3), ' Hz, \n spikes/100 ms: ', round(meanrate*0.1* length(Ecells)), sep='')
	maxbr <- max(50, ceiling(max(rate1)))
	hist(rate1, br=0: maxbr, col=IntCol, xlab='rate (Hz)', ylab='frequency', main=title)
	hist(rate1[Ecells], br=0: maxbr, add=T, col=PCCol)
	hist(rate1[actEcells], br=0: maxbr, add=T, col=aPCCol)
	
	## 2. max rate
	rate2 <- apply(ratemaps, 1, max)
	size <- rep(NA, ncells)
	for (i in 1:ncells) size[i] <- sum(ratemaps[i,,] > 0.1*rate2[i]) * dx * dx / 10000
	
	title <- paste('max rate vs. size, N =', length(actEcells))
	plot(size, rate2, pch=21, bg=IntCol, ylab='rate (Hz)', xlab='PF size (m2, r > 0.1 rmax)', main=title)
	points(size[Ecells], rate2[Ecells], pch=21, bg=PCCol)
	points(size[actEcells], rate2[actEcells], pch=21, bg=aPCCol)
	points(size[cellids], rate2[cellids], pch=21, bg=iCol)
	
	## 3. Fisher info - bound on estimation error
	errFisher <- Fisher_error(rates=ratemaps[Ecells,,], dx=dx, Tbin=0.1)
	
	K <- dim(errFisher)[1]
	L <- dim(errFisher)[2]		
	if (is.null(posmap)) pmap2 <- pmap2 <- array(1, dim=c(K, L)) else {
		pmap2 <- array(NA, dim=c(K, L))
		for (k in 1:K){
			for (l in 1:L){
				pmap2[k, l] <- (posmap[k, l] + posmap[k+1, l] + posmap[k, l+1] + posmap[k+1, l+1]) / 4
			}
		}
	}
	pmap2 <- pmap2 / sum(pmap2) 
	# mean(errFisher)
	FLB <- round(sum(pmap2 * errFisher), 2)
	title = paste('Fisher lb (100 ms), E[s] =', FLB)
	FLB <- round(sum(pmap2 * errFisher), 4)
	maxbr <- max(20, ceiling(errFisher))
	hist(errFisher, br=0:maxbr, col=PCCol, xlab='error (cm)', ylab='frequency', main=title, xlim=c(0, 25))
	
	## 4. spatial info
	info <- rep(NA, ncells)
	acells <- rep(F, ncells)
	for (i in 1:ncells) {
		if (rate1[i] > 0.1) {
			info[i] <- skaggs93.info(ratemaps[i,,], rep(1, 1600))
			acells[i] <- T
		}
	}
	info2 <- info / rate1 # bit/

	xmin <- 0.02 #min(0.1, min(info2[acells], na.rm=T))
	xmax <- max(5, ceiling(max(info2[acells], na.rm=T)))	
	ymin <- min(0.1, min(rate1[acells], na.rm=T))
	ymax <- max(50, ceiling(max(rate1[acells], na.rm=T)))	
	plot(info2[acells], rate1[acells], xlab='information (bits/sp)', ylab='firing rate (Hz)', log='xy', axes=F, pch=21, bg=IntCol, main='spatial info vs. rate', xlim=c(xmin, xmax), ylim=c(ymin,ymax))
	points(info2[Ecells], rate1[Ecells], pch=21, bg=PCCol)
	points(info2[actEcells], rate1[actEcells], pch=21, bg=aPCCol)
	axis(1, c(0.02, 0.06, 0.2, 0.6, 2, 6)); axis(2, c(0.05, 0.5, 5, 50), las=2)
	
	points(info2[cellids], rate1[cellids], pch=16, col=iCol)
	
	# plot(info, rate1, xlab='spatial information (bits/s)', ylab='firing rate (Hz)', log='xy', axes=F)
	# points(info[ddata$I_Cells], rate1[ddata$I_Cells], pch=16, col=viridis(6)[3])
	# points(info[Ecells], rate1[Ecells], pch=16, col=viridis(6,option='B')[5])
	# axis(1, c(0.02, 0.06, 0.2, 0.6, 2, 6)); axis(2, c(0.05, 0.5, 5, 50), las=2)
	# points(info[c(82, 83, 65, 86)], rate1[c(82, 83, 65, 86)], pch=16, col=viridis(6,option='B')[4])
	par(mar=c(1, 2.5, 4, 2.5))
	for (i in 1:4) image(ratemaps[cellids[i],,], col=viridis(24), main=paste('cell: ', cellids[i], ', info=', round(info2[cellids[i]], 2), ' bits/sp, rate=', round(rate1[cellids[i]], 2), 'Hz', sep=''), axes=F)
	FLB
}


estimate_posmap <- function(txy, i_data, xbrs, ybrs, sigma, graphics=F, prior=1, dim=2){
## a function to estimate the occupancy map from position data during exploration
## pos: matrix with time (sec), smoothed x and y position (cm) as columns 
## i_data: index for data used for estimation of the map
## sigma: width of a Gaussian smoothing kernel
## xbrs, ybrs: breaks for counting position 2D histogram; 
## prior: 
## return: posmap (time spent in each location, s)
	if (!dim %in% c(1,2)){
		stop('dimension should be 1 or 2')
	}
	dt <- median(diff(txy[,1]), na.rm=T)
	xw <- median(diff(xbrs))
	yw <- median(diff(ybrs))
	
	if (xw - max(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
	if (xw - min(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
	if (yw - max(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')
	if (yw - min(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')
	
	txy2 <- txy
	txy2[,2] <- txy[,2] - xbrs[1]
	ix_run <- txy2[i_data,2] %/% xw + 1

	if (dim==1){
		Nbin <- length(xbrs) - 1
		poscount <- matrix(prior, 1, Nbin)

		TT <- length(ix_run)
		for (i in 1: TT){
			poscount[1, ix_run[i]] <- poscount[1, ix_run[i]] + 1
		}
		poscount <- poscount / sum(poscount) * TT  * dt # + dt	
		pc2 <- rbind(poscount, poscount)
		pm2 <- blur(as.im(pc2), sigma/xw, normalise=T)$v
		posmap <- pm2[1,]
		if (graphics){
			plot(xbrs[-1] - xw/2, poscount, xlab='x coordinate', ylab='density', t='o', pch=15)
			lines(xbrs[-1] - xw/2, posmap, col=viridis(5)[4], lwd=2)
			points(txy[i_data,2], runif(sum(i_data), 0, max(posmap)), pch=16, cex=0.2, col=grey(0.8, alpha=0.3))
		}

	} else {
		txy2[,3] <- txy[,3] - ybrs[1]
		iy_run <- txy2[i_data,3] %/% yw + 1
	
		Nbin <- c(length(xbrs), length(ybrs)) - 1
		poscount <- matrix(prior, Nbin, Nbin)
	
		TT <- length(ix_run)
		for (i in 1: TT){
			poscount[ix_run[i], iy_run[i]] <- poscount[ix_run[i], iy_run[i]] + 1
		}
		poscount <- poscount / sum(poscount) * TT  * dt # + dt	
		posmap <- blur(as.im(poscount), sigma/xw, normalise=T)$v
	
		if (graphics){
			image(xbrs[-1] - xw/2, ybrs[-1]-yw/2, posmap, col=viridis(24), xlab='x coordinate', ylab='y coordinate')
			# image(1:20*10-5, 1:20*10-5, posmap, col=viridis(24))
			# image(1:20*10-5, 1:20*10-5, poscount, col=viridis(24))
			# print(sum(i_data))
			# print(range(i_data))
			# print(dim(txy))
			points(txy[i_data,2], txy[i_data,3], pch=16, cex=0.2, col=grey(0.8, alpha=0.3))
		}
	}	

	attr(posmap, 'n_observations') <- sum(i_data)
	attr(posmap, 'dt') <- dt
	attr(posmap, 'xcenters') <- xbrs[-1] - xw/2
	if (dim==2) attr(posmap, 'ycenters') <- ybrs[-1] - yw/2
	posmap
}


estimate_posmap_xyw <- function(xyw, xbrs, ybrs, sigma, graphics=F, prior=1, dim=2){
## a function to estimate the occupancy map from position data during exploration
## xyw: matrix with smoothed x and y position (cm) and the bin width as columns 
## xbrs, ybrs: breaks for counting position 2D histogram; 
## sigma: width of a Gaussian smoothing kernel
## prior: 
## return: posmap (time spent in each location, s)

	if (dim==1){
		xw <- median(diff(xbrs))
		
		if (xw - max(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
		if (xw - min(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
			
		Nbin <- length(xbrs) - 1
		poscount <- rep(prior, Nbin)
	
		xyw2 <- xyw
		xyw2[,1] <- xyw[,1] - xbrs[1]
		ix <- xyw2[,1] %/% xw + 1
	
		TT <- nrow(xyw)
		for (i in 1: TT){
			poscount[ix[i]] <- poscount[ix[i]] + xyw[i,3]
		}
		poscount <- poscount / sum(poscount) * sum(xyw[,3])
	
		pc2 <- rbind(poscount, poscount)
		pm2 <- blur(as.im(pc2), sigma/xw, normalise=T)$v
		posmap <- pm2[1,]

		if (graphics){
			plot(xbrs[-1] - xw/2, posmap, xlab='x coordinate', ylab='y coordinate')
		}
		attr(posmap, 'n_observations') <- nrow(xyw)
		attr(posmap, 'xcenters') <- xbrs[-1] - xw/2
	} else {
		xw <- median(diff(xbrs))
		yw <- median(diff(ybrs))
		
		if (xw - max(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
		if (xw - min(diff(xbrs)) > 1e-10) stop('xbreaks should be uniformly spaced')
		if (yw - max(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')
		if (yw - min(diff(ybrs)) > 1e-10) stop('ybreaks should be uniformly spaced')
			
		Nbin <- c(length(xbrs), length(ybrs)) - 1
		poscount <- matrix(prior, Nbin[1], Nbin[2])
	
		xyw2 <- xyw
		xyw2[,1] <- xyw[,1] - xbrs[1]
		xyw2[,2] <- xyw[,2] - ybrs[1]
	
		ix <- xyw2[,1] %/% xw + 1
		iy <- xyw2[,2] %/% yw + 1
	
		TT <- nrow(xyw)
		for (i in 1: TT){
			poscount[ix[i], iy[i]] <- poscount[ix[i], iy[i]] + xyw[i,3]
		}
		poscount <- poscount / sum(poscount) * sum(xyw[,3])
	
		posmap <- blur(as.im(poscount), sigma/xw, normalise=T)$v
	
		if (graphics){
			image(xbrs[-1] - xw/2, ybrs[-1]-yw/2, posmap, col=viridis(24), xlab='x coordinate', ylab='y coordinate')
			# image(1:20*10-5, 1:20*10-5, posmap, col=viridis(24))
			# image(1:20*10-5, 1:20*10-5, poscount, col=viridis(24))
			# print(sum(i_data))
			# print(range(i_data))
			# print(dim(txy))
			points(xyw[,1], xyw[,2], pch=16, cex=0.2, col=grey(0.8, alpha=0.3))
		}
		
		attr(posmap, 'n_observations') <- nrow(xyw)
		attr(posmap, 'xcenters') <- xbrs[-1] - xw/2
		attr(posmap, 'ycenters') <- ybrs[-1] - yw/2
	}
	posmap
}


eval_map <- function(map, xy){
	## performs bilinear interpolation to estimate the value of a map at a given point from the 4 neighbouring grid-values.
	## map: a matrix of values with attributes xcenters and ycenters with the x and y coordinates
	# map <- posmap
	# xy <- c(0.075, 0.525)

	xcenters <- attr(map, 'xcenters')
	ycenters <- attr(map, 'ycenters')
	L <- length(xcenters)
	if ((xy[1] >= min(xcenters)) & (xy[1] <= max(xcenters))){
		i1 <- max(which(xcenters <= xy[1]))
		i2 <- min(which(ycenters >= xy[1]))	
	} else {
		return(0)		
	}
	if (i1==i2){
		alpha <- 1
	} else {
		## bilinear interpolation
		dx1 <- abs(xcenters[i1] - xy[1])
		Dx <- dx1 + abs(xcenters[i2] - xy[1])
		alpha <- dx1/Dx
	}

	if ((xy[2] >= min(ycenters)) & (xy[2] <= max(ycenters))){
		j1 <- max(which(ycenters <= xy[2]))
		j2 <- min(which(ycenters >= xy[2]))
	} else {
		return(0)		
	}
	if(j1==j2) {
		beta <- 1
	} else {
		dy1 <- abs(ycenters[j1] - xy[2])
		Dy <- dy1 + abs(ycenters[j2] - xy[2])
		beta <- dy1 / Dy		
	}
		
	z1 <- (1-alpha) * map[i1,j1] + alpha * map[i2,j1]
	z2 <- (1-alpha) * map[i1,j2] + alpha * map[i2,j2]
	
	z <- (1-beta) * z1 + beta * z2
	z
}


# poscount <- matrix(0, L/xw, L/yw)
# poscount[10,10] <- 1
# posmap <- blur(as.im(poscount), 1)
# image(1:20*10-5, 1:20*10-5, posmap$v, col=viridis(24))
# plot(posmap$v[10,])
# lines(1:20, dnorm(1:20, 10, 1)/2.5, t='o')

# T_theta_min <- 1
# T_SPW_min <- 0.25

calc_poprate <- function(spt, tt, icells){
	nsp <- rep(NA, length(tt))
	dt <- median(diff(tt))
	spikes <- spt[(spt[,2] %in% icells)	,]
	N <- nrow(spikes)	

	imin <- max(which(spikes[1:1000,1] < (tt[1] - dt))) + 1
	if (imin == (-Inf)) imin <- 1

	for (i in 1:length(tt)){
		mm <- min(N, imin+1000)
		imax <- max(which(spikes[imin:mm,1] <= tt[i]))
		if (imax == (-Inf)) imax <- 0
		nsp[i] <- imax
		imin <- imin + imax
		if (imin > N) break()
	}
	# nsp <- nsp / dt
	nsp
}

get_long_periods <- function(states, N_true_min=29.5, N_false_min=7.5){
	# a function to detect periods of TRUEs and FALSEs
	# minimum duration of TRUEs: N_true_min
	# minimum interval between TRUE states: N_false_min

	# first, remove short FALSE periods by making them TRUE
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

#####################################################################

select_train_data <- function(i_run, fold=10, n=1, graphics=1){
	## select indexes of training and test data
	## we need continuous sections of run periods, so we cut at the first point when the run state changes.

	t_run <- cumsum(i_run)
	i_test <- rep(F, length(i_run))
	i_train <- rep(F, length(i_run))
	if (n == fold) {
		t_first <- round((n-1) * max(t_run) / fold)
		ii_first <- which.min((t_run - t_first)^2)
		i_cut1 <- find_cut(i_run, ii_first, 15, 'start')
		i_cut2  <- length(i_run)
		i_test[i_cut1:i_cut2] <- T
		i_train[1:(i_cut1-1)] <- T

	} else if (n == 1){
		i_cut1  <- 1
		t_last <- round((n) * max(t_run) / fold)
		ii_last <- which.min((t_run - t_last)^2)
		i_cut2 <- find_cut(i_run, ii_last, 15, 'end')
		i_test[1:i_cut2] <- T
		i_train[(i_cut2+1):length(i_run)] <- T
	} else {
		t_first <- round((n-1) * max(t_run) / fold)
		ii_first <- which.min((t_run - t_first)^2)
		i_cut1 <- find_cut(i_run, ii_first, 15, 'start')
		t_last <- round((n) * max(t_run) / fold)
		ii_last <- which.min((t_run - t_last)^2)
		i_cut2 <- find_cut(i_run, ii_last, 15, 'end')		
		i_test[i_cut1:i_cut2] <- T
		i_train[1:(i_cut1-1)] <- T
		i_train[(i_cut2+1):length(i_run)] <- T
	}


	if(graphics>0){
		layout(matrix(c(1,1,2,3), 2, byrow=T))
		tt <- 1:length(i_run)
		plot(tt, i_run, t='l', ylim=c(-.3, 1.3))
		abline(v=c(i_cut1, i_cut2), col=c(2,3), lwd=3)
	
		ii <- seq(max(1,i_cut1-30), i_cut1+30)
		plot(tt[ii], i_run[ii], t='o', ylim=c(-.3, 1.3), pch=21, bg=i_test[ii]+1)
	
		ii <- seq(i_cut2-30, min(length(tt), i_cut2+30))
		plot(tt[ii], i_run[ii], t='o', ylim=c(-.3, 1.3), pch=21, bg=i_test[ii]+1)
	}

	i_test[i_test]=i_run[i_test]
	i_train[i_train]=i_run[i_train]
	tt <- list(i_train=i_train, i_test=i_test)
	tt
		
}


find_cut <- function(i_run, ii, kmax=15, kind='start'){	
	if (kind=='start') j <- 0 else j <- 1

	k <- 1
	still_run <- TRUE
	while(still_run){
		if (i_run[ii+k] != i_run[ii]) {
			i_cut <- ii + k - j
			still_run <- FALSE
		} else if (i_run[ii-k] != i_run[ii]) {
			i_cut <- ii - k + (1 - j)
			still_run <- FALSE
		} else if (k > kmax){
			i_cut <- ii
			still_run <- FALSE
		} else {
			k <- k + 1
		}
	}
	i_cut
}


########################################################
get_theta_chunks <- function(txy, i_data, phi, t_theta, start_bins=6*pi/4, dim=2, fmin=5, fmax=12){
	# a function to calculate the theta period start and end times
	# we need continuous sections of theta during running
	# input:
	# txy: columns of time, x and y coordinates
	# i_data: index of txy with the continuous run periods indicated
	# phi: LFP theta phase
	# t_theta: time coordinates for phi
	# start_bins: theta phase for start
	#
	# output: list of theta start and end times
		
	#####################################################
	## 1. find the start time of the theta cycles
	phi <- phi %% (2*pi) - pi # transform it to the range of [-pi, pi]
	phi[which((diff(phi) < 0) & (diff(phi) > (-1 * pi)))] <- NA

	# start_bins <- 3 * pi / 4
	start_bins <- start_bins %% (2*pi) - pi # transform it to the range of [-pi, pi]
	dphi2 <- (phi - start_bins)^2
	dphi2_right_shifted <- c(NA, dphi2[1:(length(dphi2) - 1)])
	dphi2_left_shifted <- c(dphi2[2:length(dphi2)], NA)
	ii.starts <- which((dphi2 < dphi2_right_shifted) & (dphi2 < dphi2_left_shifted))
	tt.starts <- t_theta[ii.starts]
	# length(ii.starts)
	# # hist(diff(tt.starts))
	# plot(t_theta[142000:145000], phi[142000:145000])
	# abline(v=tt.starts)
	
	#####################################################
	## 2. calculate the length of the cycles - and remove extreme theta cycles, <5Hz or >12Hz
	all_theta.start.end <- cbind(tt.starts[1:(length(tt.starts)-1)], tt.starts[-1])

	binsize <- all_theta.start.end[,2] - all_theta.start.end[,1]
	shortest_bin <- 1/fmax # 1/12 corresponds to 12 Hz theta
	longest_bin <- 1/fmin # 1/5 corresponds to 5 Hz theta
	theta_bins <- (binsize > shortest_bin) & (binsize < longest_bin)
	true_theta.start.end <- all_theta.start.end[theta_bins,]
	
	#####################################################
	## 3. position indices during theta cycling
	run_theta_OK <- approx(all_theta.start.end[,1], theta_bins, txy[,1], method='constant', f=0)$y # 1 if theta is OK
	run_theta_OK[is.na(run_theta_OK)] <- 0 # outside the theta intervals...
	i_data[!run_theta_OK] <- FALSE
	i_data_long <- get_long_periods(i_data, N_true_min=30, N_false_min=0)

	# hist(true_theta.start.end[,2] - true_theta.start.end[,1])
	
	#####################################################
	## 4. find theta during running
	## the start of the theta is during run
	theta_run1 <- approx(txy[,1], i_data_long, true_theta.start.end[,1], method='constant', f=0)$y
	theta.runstart <- true_theta.start.end[!!theta_run1,]
	## the end of the theta is during run
	theta_run2 <- approx(txy[,1], i_data_long, theta.runstart[,2], method='constant', f=1)$y
	theta.start.end <- theta.runstart[!!theta_run2,]


	#########################################################
	## 5. find the chunks - start and end of theta&run periods, in run-coordinates
	i_chunks_start_end <- find_chunks(i_data_long)
	# print(i_chunks_start_end[1,])
	
	if (is.null(dim(i_chunks_start_end))){
		t_chunks_start_end <- matrix(txy[i_chunks_start_end], 1)
		i_chunks_start_end <- matrix(i_chunks_start_end, 1)
	} else {
		t_chunks_start_end <- i_chunks_start_end
		t_chunks_start_end[,1] <- txy[i_chunks_start_end[,1],1]
		t_chunks_start_end[,2] <- txy[i_chunks_start_end[,2],1]		
	}
	# print(t_chunks_start_end[1,])
	# abline(v=t_chunks_start_end[1,], col=c(2,3), lwd=2)		
	
	theta <- list() # list with the theta chunks
	n.chunks <- nrow(i_chunks_start_end)
	i.chunk <- 1
		
	for (i.chunk in 1:n.chunks){
		tmin <- t_chunks_start_end[i.chunk,1]
		tmax <- t_chunks_start_end[i.chunk,2]
		i_theta_chunks <- which((theta.start.end[,1] > tmin) & (theta.start.end[,1] < tmax))
		theta_chunks <- theta.start.end[i_theta_chunks,]
		theta[[i.chunk]] <- theta_chunks
	}
	# print(theta[[1]])
	theta
}


########################################################
get_phase_chunks <- function(sp, txy, i_data, phi=NULL, t_theta=NULL, start_phi=NULL, thetas = NULL, L_phi, cellids, dim=2, fmin=5, fmax=12){
	# a function to create spike-count chunks for decoding during theta
	# spikes are counted in windows proportional to the duration of the theta cycle
	#
	# input:
	# sp: two columns: spikes, cell
	# txy: columns of time, x and y coordinates
	# i_data: index of txy with the continuous run periods indicated
	# phi: LFP theta phase
	# t_theta: time coordinates for phi
	# start_phi: theta phase for start
	# L_phi: phase length L
	# cellids: cells to be included in the rasters
	#
	# output: list of phase chunks
	# each chunk being a matrix of spike counts for all cells in cellids
	# the position of the animal is also provided for each chunk for the same time bins
	# sp <- ddata$spikes;  txy <- txyv; i_data <- i_run; phi <- h_theta; t_theta <- t_theta; start_phi <- 6*pi/4; L_phi=2*pi/3; cellids <- act_Ecells

	#########################################################
	## first, find the theta chunks - start and end of theta periods	and cycles during run
	#########################################################
	#########################################################
	# find the required theta phases
	if (is.null(thetas)){
		if (is.null(phi)) stop('error: theta phase must be provided!')
		if (is.null(t_theta)) stop('error: theta time must be provided!')
		if (is.null(start_phi)) stop('error: theta start phase must be provided!')
		thetas <- get_theta_chunks(txy, i_data, phi, t_theta, start_bins=start_phi, dim=dim, fmin=fmin, fmax=fmax)
	}

	################################################################
	
	chunks <- list()
	n.chunks <- length(thetas)
	colnames(sp) <- c('time', 'cell')
	N.cells <- max(sp[,2])
	ii.cells <- sp[,2] %in% cellids
	ii.excluded <- which(!(1:N.cells %in% cellids))
	sp <- sp[ii.cells,]
	
	for (i.chunk in 1:n.chunks){
		theta_chunks <- thetas[[i.chunk]]
		tmin <- min(theta_chunks) - 0.1
		tmax <- max(theta_chunks) + 0.1

		bin_width <- (theta_chunks[,2] - theta_chunks[,1]) * L_phi / (2*pi)
		theta_ends <- theta_chunks[,1] + bin_width
		theta_chunks[,2] <- theta_ends
		
		tt <- rowMeans(theta_chunks)
		nbins <- length(tt)
		xpos <- approx(txy[,1], txy[,2], xout=tt)$y
		if (dim == 2) {
			ypos <- approx(txy[,1], txy[,3], xout=tt)$y
			dir <- approx(txy[,1], txy[,5], xout=tt)$y
		}
		
		i_sp <- which((sp[,1] > tmin) & (sp[,1]<tmax))
		if (length(i_sp) > 0){
			sp_i <- sp[i_sp,]
	
			rast <- spt2raster_bins(sp_i, theta_chunks, N.cells)
			rast[,ii.excluded] <- NA
			
			if (dim==2){
				chunks[[i.chunk]] <- list(rast=rast, pos=cbind(tt,xpos, ypos, dir, bin_width))
			} else {
				chunks[[i.chunk]] <- list(rast=rast, pos=cbind(tt,xpos, bin_width))				
			}
			# plot(xpos, ypos, xlim=c(0, 200), ylim=c(0, 200))
			# readline(nrow(sp_i))
			# cat(i.chunk, ' ')			
		} else {
			if (dim==2){
				chunks[[i.chunk]] <- list(rast=NA, pos=cbind(tt,xpos, ypos, dir, bin_width))
			} else {
				chunks[[i.chunk]] <- list(rast=NA, pos=cbind(tt, xpos, bin_width))				
			}
		}
	}	
	chunks
}


########################################################
get_20ms_seqs <- function(sp, txy, i_data, phi=NULL, t_theta=NULL, start_bins=NULL, cellids, thetas = NULL, dim=2, fmin=5, fmax=12){
	# a function to create spike-count chunks for decoding during theta
	# each chunk is divided into 20ms bins parts and separate chunks are provided for each with similar spike counts
	#
	# input:
	# sp: two columns: spikes, cell
	# txy: columns of time, x and y coordinates
	# i_data: index of txy with the continuous run periods indicated
	# phi: LFP theta phase
	# t_theta: time coordinates for phi
	# start_bins: theta phase for start
	# cellids: cells to be included in the rasters
	# optional:
	# thetas: the time 
	#
	# output: list of theta chunks
	# each chunk being a list with 
	#		- matrix of spike counts in 20 ms overlapping bins for each cells
	# 		- the bin start, end, theta czcle, position and direction
	# sp <- spt;  txy <- txyvd; i_data <- i_run; phi <- h_theta; start_phi <- theta_start; cellids <- act_Ecells


	#########################################################
	# find the required theta phases
	if (is.null(thetas)){
		if (is.null(phi)) stop('error: theta phase must be provided!')
		if (is.null(t_theta)) stop('error: theta time must be provided!')
		if (is.null(start_bins)) stop('error: theta start phase must be provided!')
		thetas <- get_theta_chunks(txy, i_data, phi, t_theta, start_bins=start_bins, dim=dim, fmin=fmin, fmax=fmax)
	}

	################################################################
	
	rasts <- list()
	n.chunks <- length(thetas)
	N.cells <- max(sp[,2])
	ii.cells <- sp[,2] %in% cellids
	sp <- sp[ii.cells,]
	colnames(sp) <- c('time', 'cell')
	ii.excluded <- which(!(1:N.cells %in% cellids))
	
	for (i.chunk in 1:n.chunks){
		theta_chunks <- thetas[[i.chunk]]
		n.cycles <- nrow(theta_chunks)

		tmin <- min(theta_chunks) 
		nbins <- floor(diff(range(theta_chunks)) / 0.005) # start at tmin and use 5 ms bins 
		tmax <- tmin + (nbins+1) * 0.005
		
		bin_start <- seq(tmin, by=0.005, length=nbins)
		bin_end <- bin_start + 0.02
		i.cycle <- rep(0, length(bin_start))
		for (ii in 1: n.cycles) {
			i.cycle[bin_start >= theta_chunks[ii,1]] <- ii
		}
		bin_start_end <- cbind(bin_start, bin_end, i.cycle)
		
		ii_sp <- (sp[,1] >= tmin) & (sp[,1] <= tmax)
		sp_i <- sp[ii_sp,]
		
		rast <- spt2raster_overlapping_bins(sp_i, bin_start_end, N.cells)
		if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
		bin_mids <- rowMeans(bin_start_end[,1:2])				

		t_start <- c(theta_chunks[,1], max(theta_chunks[,2]))
		t_phi <- seq(0, by=2*pi, length=length(t_start))
		phi_bins <- approx(t_start, t_phi, bin_start)$y 

		xpos <- approx(txy[,1], txy[,2], xout=bin_mids)$y
		if (dim == 2) {
			ypos <- approx(txy[,1], txy[,3], xout=bin_mids)$y
			dir <- approx(txy[,1], txy[,5], xout=bin_mids)$y
			rasts[[i.chunk]] <- list(rast=rast, pos=cbind(bin_start_end, xpos, ypos, dir, phi_bins))
		} else {
			rasts[[i.chunk]] <- list(rast=rast, pos=cbind(bin_start_end, xpos))			
		}
	cat('chunk: ', i.chunk, '\n')
	}	
	rasts
}

########################################################
get_theta_seqs <- function(sp, txy, i_data, phi=NULL, t_theta=NULL, start_bins=NULL, cellids, thetas = NULL, dim=2, fmin=5, fmax=12, equal_spikes=T){
	# a function to create spike-count chunks for decoding during theta
	# each theta cycle is divided into 3 parts and separate chunks are provided for each with similar spike counts
	#
	# input:
	# sp: two columns: spikes, cell
	# txy: columns of time, x and y coordinates
	# i_data: index of txy with the continuous run periods indicated
	# phi: LFP theta phase
	# t_theta: time coordinates for phi
	# start_bins: theta phase for start
	# cellids: cells to be included in the rasters
	# optional:
	# thetas: list, 
	# 		- each element is a matrix corresponding to continuous theta periods
	# 		- each row contains the start and end of a theta cycle
	# output: list of theta chunks
	# each chunk being a matrix of spike counts in each theta cycle for all cells in cellids
	# the position of the animal is also provided for each chunk in the same time bins
	# note, that the time resolution is not constant, it changes with the theta periods
	# sp <- ddata$spikes;  txy <- txyv; i_data <- i_run; phi <- h_theta; t_theta <- t_theta; start_bins <- 6*pi/4; cellids <- act_Ecells

	#########################################################
	# find the required theta phases
	if (is.null(thetas)){
		if (is.null(phi)) stop('error: theta phase must be provided!')
		if (is.null(t_theta)) stop('error: theta time must be provided!')
		if (is.null(start_bins)) stop('error: theta start phase must be provided!')
		thetas <- get_theta_chunks(txy, i_data, phi, t_theta, start_bins=start_bins, dim=dim, fmin=fmin, fmax=fmax)
	}

	################################################################
	
	rast_1 <- list()
	rast_2 <- list()
	rast_3 <- list()

	n.chunks <- length(thetas)
	N.cells <- max(sp[,2])
	ii.cells <- sp[,2] %in% cellids
	sp <- sp[ii.cells,]
	colnames(sp) <- c('time', 'cell')
	ii.excluded <- which(!(1:N.cells %in% cellids))
	
	for (i.chunk in 1:n.chunks){
		theta_chunks <- thetas[[i.chunk]]
		tmin <- min(theta_chunks) - 0.1
		tmax <- max(theta_chunks) + 0.1

		i_sp <- which((sp[,1] >= tmin) & (sp[,1] <= tmax))
		spikes <- sp[i_sp,]
		# spikes <- spikes_i[which(!(spikes_i[,2] %in% ii.excluded)),]

		if (length(spikes) > 0){
			if (equal_spikes){
				rast123 <- sp2thetaRasts(spikes, theta_chunks, N_bins=3, N.cells)
			} else {
				rast123 <- list()
				bin_size <- (theta_chunks[,2] - theta_chunks[,1]) / 3
				bin_start_end <- cbind(theta_chunks[,1], theta_chunks[,1] + bin_size)
				bin_mids <- rowMeans(bin_start_end)				
				rast123$raster1 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123$bintime1 <- cbind(bin_size, bin_start_end)
				
				bin_start_end <- cbind(theta_chunks[,1] + bin_size, theta_chunks[,1] + 2*bin_size)
				bin_mids <- rowMeans(bin_start_end)				
				rast123$raster2 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123$bintime2 <- cbind(bin_size, bin_start_end)

				bin_start_end <- cbind(theta_chunks[,1] + 2* bin_size, theta_chunks[,2])
				bin_mids <- rowMeans(bin_start_end)				
				rast123$raster3 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123$bintime3 <- cbind(bin_size, bin_start_end)
			}

			rast <- rast123$raster1
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123$bintime1[,2])$y
			bin_width <- rast123$bintime1[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123$bintime1[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123$bintime1[,2])$y		
				rast_1[[i.chunk]] <- list(rast=rast, pos=cbind(rast123$bintime1[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_1[[i.chunk]] <- list(rast=rast, pos=cbind(rast123$bintime1[,2], xpos, bin_width))			
			}

			rast <- rast123$raster2
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123$bintime2[,2])$y
			bin_width <- rast123$bintime2[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123$bintime2[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123$bintime2[,2])$y		
				rast_2[[i.chunk]] <- list(rast=rast, pos=cbind(rast123$bintime2[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_2[[i.chunk]] <- list(rast=rast, pos=cbind(rast123$bintime2[,2], xpos, bin_width))				
			}

			rast <- rast123$raster3
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123$bintime3[,2])$y
			bin_width <- rast123$bintime3[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123$bintime3[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123$bintime3[,2])$y		
				rast_3[[i.chunk]] <- list(rast=rast, pos=cbind(rast123$bintime3[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_3[[i.chunk]] <- list(rast=rast, pos=cbind(rast123$bintime3[,2], xpos, bin_width))				
			}

		} else {
			if (dim == 2) {
				rast_1[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))			
				rast_2[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))			
				rast_3[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))
			} else {
				rast_1[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))			
				rast_2[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))			
				rast_3[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))				
			}
		}
	}	
	
	
	rasts <- list(past=rast_1, present=rast_2, future=rast_3)
	rasts
}


########################################################
get_theta_seqs6 <- function(sp, txy, i_data, phi=NULL, t_theta=NULL, start_bins=NULL, cellids, thetas = NULL, dim=2, fmin=5, fmax=12, equal_spikes=T){
	# a function to create spike-count chunks for decoding during theta
	# each theta cycle is divided into 3 parts and separate chunks are provided for each with similar spike counts
	#
	# input:
	# sp: two columns: spikes, cell
	# txy: columns of time, x and y coordinates
	# i_data: index of txy with the continuous run periods indicated
	# phi: LFP theta phase
	# t_theta: time coordinates for phi
	# start_bins: theta phase for start
	# cellids: cells to be included in the rasters
	# optional:
	# thetas: the time 
	#
	# output: list of theta chunks
	# each chunk being a matrix of spike counts in each theta cycle for all cells in cellids
	# the position of the animal is also provided for each chunk in the same time bins
	# note, that the time resolution is not constant, it changes with the theta periods
	# sp <- ddata$spikes;  txy <- txyv; i_data <- i_run; phi <- h_theta; t_theta <- t_theta; start_bins <- 6*pi/4; cellids <- act_Ecells

	#########################################################
	# find the required theta phases
	if (is.null(thetas)){
		if (is.null(phi)) stop('error: theta phase must be provided!')
		if (is.null(t_theta)) stop('error: theta time must be provided!')
		if (is.null(start_bins)) stop('error: theta start phase must be provided!')
		thetas <- get_theta_chunks(txy, i_data, phi, t_theta, start_bins=start_bins, dim=dim, fmin=fmin, fmax=fmax)
	}

	################################################################
	
	rast_1 <- list()
	rast_2 <- list()
	rast_3 <- list()
	rast_4 <- list()
	rast_5 <- list()
	rast_6 <- list()

	n.chunks <- length(thetas)
	N.cells <- max(sp[,2])
	ii.cells <- sp[,2] %in% cellids
	sp <- sp[ii.cells,]
	colnames(sp) <- c('time', 'cell')
	ii.excluded <- which(!(1:N.cells %in% cellids))
	
	for (i.chunk in 1:n.chunks){
		theta_chunks <- thetas[[i.chunk]]
		tmin <- min(theta_chunks) - 0.1
		tmax <- max(theta_chunks) + 0.1

		i_sp <- which((sp[,1] >= tmin) & (sp[,1] <= tmax))
		spikes <- sp[i_sp,]
		# spikes <- spikes_i[which(!(spikes_i[,2] %in% ii.excluded)),]

		if (length(spikes) > 0){
			if (equal_spikes){
				rast123456 <- sp2thetaRasts(spikes, theta_chunks, N_bins=6, N.cells)
			} else {
				rast123456 <- list()
				bin_size <- (theta_chunks[,2] - theta_chunks[,1]) / 6
				bin_start_end <- cbind(theta_chunks[,1], theta_chunks[,1] + bin_size)
				bin_mids <- rowMeans(bin_start_end)				
				rast123456$raster1 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123456$bintime1 <- cbind(bin_size, bin_start_end)
				
				bin_start_end <- cbind(theta_chunks[,1] + bin_size, theta_chunks[,1] + 2*bin_size)
				bin_mids <- rowMeans(bin_start_end)				
				rast123456$raster2 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123456$bintime2 <- cbind(bin_size, bin_start_end)

				bin_start_end <- cbind(theta_chunks[,1] + 2*bin_size, theta_chunks[,1] + 3*bin_size)
				bin_mids <- rowMeans(bin_start_end)				
				rast123456$raster3 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123456$bintime3 <- cbind(bin_size, bin_start_end)

				bin_start_end <- cbind(theta_chunks[,1] + 3*bin_size, theta_chunks[,1] + 4*bin_size)
				bin_mids <- rowMeans(bin_start_end)				
				rast123456$raster4 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123456$bintime4 <- cbind(bin_size, bin_start_end)

				bin_start_end <- cbind(theta_chunks[,1] + 4*bin_size, theta_chunks[,1] + 5*bin_size)
				bin_mids <- rowMeans(bin_start_end)				
				rast123456$raster5 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123456$bintime5 <- cbind(bin_size, bin_start_end)

				bin_start_end <- cbind(theta_chunks[,1] + 5* bin_size, theta_chunks[,2])
				bin_mids <- rowMeans(bin_start_end)				
				rast123456$raster6 <- spt2raster_bins(spikes, bin_start_end, N.cells)
				rast123456$bintime6 <- cbind(bin_size, bin_start_end)
			}

			rast <- rast123456$raster1
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123456$bintime1[,2])$y
			bin_width <- rast123456$bintime1[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123456$bintime1[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123456$bintime1[,2])$y		
				rast_1[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime1[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_1[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime1[,2], xpos, bin_width))			
			}

			rast <- rast123456$raster2
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123456$bintime2[,2])$y
			bin_width <- rast123456$bintime2[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123456$bintime2[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123456$bintime2[,2])$y		
				rast_2[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime2[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_2[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime2[,2], xpos, bin_width))				
			}

			rast <- rast123456$raster3
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123456$bintime3[,2])$y
			bin_width <- rast123456$bintime3[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123456$bintime3[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123456$bintime3[,2])$y		
				rast_3[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime3[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_3[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime3[,2], xpos, bin_width))				
			}


			rast <- rast123456$raster4
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123456$bintime4[,2])$y
			bin_width <- rast123456$bintime4[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123456$bintime4[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123456$bintime4[,2])$y		
				rast_4[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime4[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_4[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime4[,2], xpos, bin_width))				
			}


			rast <- rast123456$raster5
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123456$bintime5[,2])$y
			bin_width <- rast123456$bintime5[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123456$bintime5[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123456$bintime5[,2])$y		
				rast_5[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime5[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_5[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime5[,2], xpos, bin_width))				
			}


			rast <- rast123456$raster6
			if (length(ii.excluded) > 0) rast[,ii.excluded] <- NA
			xpos <- approx(txy[,1], txy[,2], xout=rast123456$bintime6[,2])$y
			bin_width <- rast123456$bintime6[,1]
			if (dim == 2) {
				ypos <- approx(txy[,1], txy[,3], xout=rast123456$bintime6[,2])$y
				dir <- approx(txy[,1], txy[,5], xout=rast123456$bintime6[,2])$y		
				rast_6[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime6[,2], xpos, ypos, dir, bin_width))
			} else {
				rast_6[[i.chunk]] <- list(rast=rast, pos=cbind(rast123456$bintime6[,2], xpos, bin_width))				
			}

		} else {
			if (dim == 2) {
				rast_1[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))			
				rast_2[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))			
				rast_3[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))
				rast_4[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))
				rast_5[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))
				rast_6[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA, NA, NA))
			} else {
				rast_1[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))			
				rast_2[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))			
				rast_3[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))				
				rast_4[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))			
				rast_5[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))			
				rast_6[[i.chunk]] <- list(rast=NA, pos=cbind(NA, NA, NA))				
			}
		}
	}	
	
	
	rasts <- list(r1=rast_1, r2=rast_2, r3=rast_3, r4=rast_4, r5=rast_5, r6=rast_6)
	rasts
}


########################################################
## a combined plot to visualise the theta sequences
########################################################
plot_theta_seqs <- function(spt, theta_i, txy, ppf, ratemaps, cellids){
	# spt: spike times and cell ids
	# thetas: theta period boundaries for the given chunk
	# txy: position data for the given chunk
	# ppf: decoded posterior mean past, present and future for the given chunk
	# ratemaps: smoothed ratemaps for all cells
	# output: 
	# plot of the trajectory at this chunk, the current position and the decoded segments, the time and ID of the spikes and the ratemaps
	# 1. trajectory
	dev.set(2)
	par(mar=c(1,1,1,1))
	N <- nrow(theta_i)
	ttime <- rowMeans(theta_i)
	xpos <- approx(txy[,1], txy[,2], xout=ttime)$y
	ypos <- approx(txy[,1], txy[,3], xout=ttime)$y

	#########################################################
	# 2. spikes
	colnames(spt) <- c('time', 'cell')
	N.cells <- max(spt[,2])
	ii.cells <- spt[,2] %in% cellids
	spt <- spt[ii.cells,]
	
	tmin <- min(theta_i) - 0.1
	tmax <- max(theta_i) + 0.1
	ii.spikes <- which((spt[,1] > tmin) & (spt[,1] < tmax))
	spikes <- spt[ii.spikes,]

	for (i in 1:N){
	
		# I. plot trajectory
		dev.set(2)
		plot(xpos, ypos, axes=F, xlim=range(0, 200), ylim=range(0,200), pch=1, col=viridis(N), xlab='', ylab='')
		points(xpos[i], ypos[i], pch=16, col=viridis(N)[i])
		points(ppf[i,1,], ppf[i,2,], t='o', lwd=1, col=viridis(N, option='D')[i], pch=16:18, cex=spike_cex[i])
		box()
		# axis(1)
		# axis(2)
		
		# 	II. plot spikes
		rates_i <- ratemaps[,xpos[i] %/% 5, ypos[i] %/% 5]
		cellsort <- sort(sort(rates_i, index=T)$ix, index=T)$ix
		
		if (i ==1) t0 <- theta_i[i,1] - 0.1 else t0 <- theta_i[i-1,1]
		t1 <- theta_i[i,1]
		t2 <- theta_i[i,2]
		if (i==N) t3 <- theta_i[i,2] + 0.1 else t3 <- theta_i[i+1,2]
	
		iii.spikes <- which((spikes[,1] > t0) & (spikes[,1] < t3))
		spikes_i <- spikes[iii.spikes,]
		base <- 11
		cols <- colormap(colormap=colormaps$jet, N.cells)[seq(1, by=base, length=N.cells) %% N.cells]
		dev.set(3)
		plot(spikes_i[,1], cellsort[spikes_i[,2]], col=cols[cellsort[spikes_i[,2]]], pch=18, axes=F, xlab='time (s)', ylab='cells', ylim=c(0, N.cells))
		axis(1); axis(2, las=2)
		abline(v=c(t1, t2))
		
		# III. plot ratemaps
		N.active <- length(unique(spikes_i[,2]))
		kk <- sort(unique(spikes_i[,2]))
		map_order <- sort(rates_i[kk], index=T)$ix		
		
		Ncol <- 6
		Nrow <- (N.active - 1) %/% 6 + 1
		xx <- seq(2.5, 197.5, by=5)
		dev.set(4)
		par(mfcol = c(Nrow, Ncol))
		par(mar=c(1,1,1,1))
		for (k in map_order){
			title <- paste(cellsort[kk[k]], ', ', sum(spikes_i[,2] == kk[k]), ', ', round(max(ratemaps[kk[k],,])), sep='')
			image(xx, xx, ratemaps[kk[k],,], col=viridis(36), axes=F, main='')
			mtext(title, 3, col=cols[cellsort[kk[k]]], font=2, cex=0.8)
			points(xpos[i], ypos[i], pch=16, col=rgb(1, 0.5, 0))
			
		}
		readline(i)
	}
	
}

########################################################
run_chunks <- function(sp, txy, i_data, dt, cellids, tshift=dt){
	# a function to create spike-count chunks for decoding during theta
	# output: list of chunks
	# each chunk being a matrix of spike counts in dt bins for all cells in cellids
	# the position of the animal is also provided for each chunk with the same temporal resolution	
	# tshift: time difference between neighbouring bins for spike counts
	# i_data: is a binary vector indexing the time bins to be included in the analysis. 

	## first, find the chunks - start and end of run periods	
	i_chunks_start_end <- find_chunks(i_data)
	
	if (is.null(dim(i_chunks_start_end))){
		t_chunks_start_end <- matrix(txy[i_chunks_start_end], 1)
		i_chunks_start_end <- matrix(i_chunks_start_end, 1)
	} else {
		t_chunks_start_end <- i_chunks_start_end
		t_chunks_start_end[,1] <- txy[i_chunks_start_end[,1],1]
		t_chunks_start_end[,2] <- txy[i_chunks_start_end[,2],1]		
	}
	
	# plot(i_data[1:5000], t='l')
	# abline(v=i_chunks_start_end[,1], col=3)
	# abline(v=i_chunks_start_end[,2], col=2)

	if ((dt / tshift) != round(dt/tshift)){
		tshift <- dt / round(dt/tshift)
		print(paste('dt is not an integer multiple of tshift. We change tshift to '), tshift)		
	}
	
	chunks <- list()
	n.chunks <- nrow(i_chunks_start_end)
	colnames(sp) <- c('time', 'cell')
	N.cells <- max(sp[,2])
	ii.cells <- sp[,2] %in% cellids
	ii.excluded <- which(!(1:N.cells %in% cellids))
	sp <- sp[ii.cells,]
	
	ncells_merge <- dt/tshift ## 
	if (ncells_merge%%2==0) even <- T else even <- F

	ii.chunk <- 1
	for (i.chunk in 1:n.chunks){
		tmin <- floor(t_chunks_start_end[i.chunk,1] / tshift) * tshift
		tmax <- ceiling(t_chunks_start_end[i.chunk,2] / tshift) * tshift
		if ((tmax - tmin) > (2*dt)){
			if (even){
				tt_pos <- seq.int(tmin/tshift, tmax/tshift)*tshift # when we merge 2 cells, the boundary between the bins is the right time coordinate
				Tmin <- tmin - tshift*ncells_merge/2 # Tmin and Tmax for the spikes - we calculate bins by tshift and merge them accordingly 
				Tmax <- tmax + tshift*ncells_merge/2
			} else {
				tt_pos <- seq.int(tmin/ tshift, tmax / tshift-1) * tshift + tshift/2	 # if we merge odd (1,3,5) cells, then their center is the relevant
				Tmin <- tmin - tshift*(ncells_merge-1)/2 # Tmin and Tmax for the spikes - we calculate bins by tshift and merge them accordingly 
				Tmax <- tmax + tshift*(ncells_merge-1)/2
			}
			
			xpos <- approx(txy[,1], txy[,2], xout= tt_pos, rule=2)$y
			ypos <- approx(txy[,1], txy[,3], xout= tt_pos, rule=2)$y
			
			i_sp <- which((sp[,1] > Tmin) & (sp[,1]<Tmax))
			sp_i <- sp[i_sp,]
			if (nrow(sp_i) > 1){
				rast <- spt2raster(sp_i, dt, Tmin, Tmax, N.cells, tshift)
				rast[,ii.excluded] <- NA
				
			
				chunks[[ii.chunk]] <- list(rast=rast, pos=cbind(xpos, ypos))
				ii.chunk <- ii.chunk + 1
				# plot(xpos, ypos, xlim=c(0, 200), ylim=c(0, 200))
				# readline(nrow(sp_i))
				if (length(xpos) != nrow(rast)) warning(paste('raster and position is not matched! position:', length(xpos), ', raster:', nrow(rast)))
				cat(i.chunk, ' ')
			} 
		} else {
			cat(' s ')
		}
	}	
	chunks
}


find_chunks <- function(i_data){
	ii_start <- i_data
	ii_end <- c(i_data, F)
	
	i_start <- min(which(ii_start == T))
	ii_end[1:i_start] <- T	
	i_end <- min(which(ii_end == F)) - 1
	ii_start[1:i_end] <- F
	i_chunks_start_end <- c(i_start, i_end)

	while(sum(ii_start) > 0){
		i_start <- min(which(ii_start == T))
		ii_end[1:i_start] <- T	
		i_end <- min(which(ii_end == F)) - 1
		ii_start[1:i_end] <- F
		i_chunks_start_end <- rbind(i_chunks_start_end, c(i_start, i_end))
	}

	i_chunks_start_end
}
	
## convert spike times to rasters for decoding
## by default we use uniformly spaced time bins between Tmin and Tmax, with size dt and shift tshift
spt2raster <- function(spt, dt, Tmin=NULL, Tmax=NULL, N.cells=NULL, tshift=dt){
	if (!(colnames(spt)[1] %in% c("time", "cell"))) stop("columns must be named as time and cell")
	if (!(colnames(spt)[2] %in% c("time", "cell"))) stop("columns must be named as time and cell")
	
	if ((dt / tshift) != round(dt/tshift)){
		tshift <- dt / round(dt/tshift)
		print(paste('dt is not an integer multiple of tshift. We change tshift to '), tshift)		
	}	
	
	if (is.null(N.cells)) {
		N.cells <- length(unique(spt[,"cell"]))
		ii.cells <- unique(spt[,"cell"])
		renumber_cells <- T
	} else renumber_cells <- F

	if (is.null(Tmin)) {
		Tmin <- min(spt[,"time"])
	} else {
		spt <- spt[spt[,"time"]>Tmin,]		
	}
	if (is.null(Tmax)) {
		Tmax <- max(spt[,"time"]) 
	} else {
		spt <- spt[spt[,"time"]<Tmax,]
	}
	spt[,"time"] <- spt[,"time"] - Tmin
	spt[spt[,"time"]==0,"time"] <- tshift / 2

	Tmax <- Tmax - Tmin
	
	Lmax <- ceiling((Tmax / tshift))
	NN <- Lmax * N.cells
	N.spikes <- nrow(spt)
	sparseness <- N.spikes / NN
	if (sparseness < 1/50)	rast <- Matrix(0, Lmax, N.cells)	else rast <- matrix(0, Lmax, N.cells)
	for (i in 1:N.spikes){
		i.sp <- ceiling(spt[i,"time"] / tshift)
		if (renumber_cells) i.cell <- which( ii.cells == spt[i,"cell"]) else i.cell <- spt[i, "cell"]
		rast[i.sp, i.cell] <- rast[i.sp, i.cell] + 1
	}
	
	if (tshift != dt){
		ncells_merge <- dt/tshift # number of cells to merge
		n_side <- ncells_merge - 1 # we loose this many cells because of merging...
		LLmax <- Lmax - n_side
		if (sparseness < 1/100)	rast2 <- Matrix(0, LLmax, N.cells)	else rast2 <- matrix(0, LLmax, N.cells)	
		for (i in 1:LLmax){
			rast2[i,] <- colSums(rast[i:(i+ ncells_merge-1),])
		}
		rast <- rast2
	}
	rast
}

## convert spike times to rasters for decoding
## we use timebins provided 
spt2raster_bins <- function(spt, bin.start.end, N.cells=NULL){
	## put spikes in bis defined in bin.start.end
	## spt: two columns, with the name 'time' (ordered) and 'cell'.
	## bin.start.end has 2 columns, 
	##		1. start of the bins
	## 	2. end of the bins
	## time in both collumns should be increasing
	## bins can NOT overlap

	if (!(colnames(spt)[1] %in% c("time", "cell"))) stop("columns must be named as time and cell")
	if (!(colnames(spt)[2] %in% c("time", "cell"))) stop("columns must be named as time and cell")
	if (min(diff(bin.start.end[,1])) < 0) stop('bin starts should be increasing')
	if (min(diff(bin.start.end[,2])) < 0) stop('bin ends should be increasing')

	nbins <- nrow(bin.start.end)
	if (sum((bin.start.end[2:nbins,1] - bin.start.end[1:(nbins-1),2]) < 0) > 0) stop('time bins are ovelapping')
	if (max(diff(spt[,'time'])) < 0) {
		print('spike times should be increasing')
		ix <- sort(spt[,'time'], index.ret=T)$ix
		spt <- spt[ix,]
	}
		
	if (is.null(N.cells)) {
		N.cells <- length(unique(spt[,"cell"]))
		ii.cells <- unique(spt[,"cell"])
		renumber_cells <- T
	} else renumber_cells <- F
	
	NN <- nbins * N.cells

	first_spike <- min(which(spt[,"time"] > bin.start.end[1,1]))
	last_spike <- max(which(spt[,"time"] < bin.start.end[nbins,2]))
	spt_pos <- spt[first_spike:last_spike,, drop=F]

	N.spikes <- nrow(spt_pos)
	
	# sparseness <- N.spikes * sum(bin.start.end[,2] - bin.start.end[,1]) / (max(bin.start.end) - min(bin.start.end)) / NN
	# if (sparseness < 1/50)	rast <- Matrix(0, nbins, N.cells)	else rast <- matrix(0, nbins, N.cells)
	rast <- matrix(0, nbins, N.cells)
	
	i_start <- 1	
	i_end <- min(nbins, i_start + 10)
	for (i in 1:N.spikes){
		t.sp <- spt_pos[i,"time"]
		
		bin_not_found <- T
		while(bin_not_found){
			i.bin <- i_start -1 + max(which(bin.start.end[i_start:i_end,1] <= t.sp))
			if (i.bin > 0) bin_not_found <- F else i_end <- min(nbins, i_end + 10)
		}
		if (t.sp < bin.start.end[i.bin,2]){
			if (renumber_cells) i.cell <- which( ii.cells == spt_pos[i,"cell"]) else i.cell <- spt_pos[i, "cell"]
			rast[i.bin, i.cell] <- rast[i.bin, i.cell] + 1
		}
		i_start <- i.bin
		i_end <- min(nbins, i_start + 10)
	}
	rast
}

## function for overlapping bins
spt2raster_overlapping_bins <- function(spikes, bin.start.end, N.cells=NULL){
	## put spikes in bis defined in bin.start.end
	## spikes: two columns, with the name 'time' (ordered) and 'cell'.
	## bin.start.end has 2 columns, 
	##		1. start of the bins
	## 	2. end of the bins
	## time in both collumns should be increasing
	## all timepoints in the range(bin.start.end) should be in at least one of the bins
	## bins can OVERLAP
	
	if (!(colnames(spikes)[1] %in% c("time", "cell"))) stop("columns must be named as time and cell")
	if (!(colnames(spikes)[2] %in% c("time", "cell"))) stop("columns must be named as time and cell")

	if (min(diff(bin.start.end[,1])) < 0) stop('bin starts should be increasing')
	if (min(diff(bin.start.end[,2])) < 0) stop('bin ends should be increasing')
	# if (sum((bin.start.end[2:nbins,1] - bin.start.end[1:(nbins-1),2]) < 0) > 0) ovelapping <- T

	nbins <- nrow(bin.start.end)
	if (max(bin.start.end[2:nbins,1] - bin.start.end[1:(nbins-1),2]) > 1e-15) stop('bins are not tiling the whole interval!')
	if (max(diff(spikes[,'time'])) < 0) {
		print('spike times should be increasing')
		ix <- sort(spikes[,'time'], index.ret=T)$ix
		spikes <- spikes[ix,]
	}
	
	if (is.null(N.cells)) {
		N.cells <- length(unique(spikes[,"cell"]))
		ii.cells <- unique(spikes[,"cell"])
		renumber_cells <- T
	} else renumber_cells <- F
	
	bin_width <- median(bin.start.end[,2] - bin.start.end[,1])
	bin_change <- median(diff(bin.start.end[,1]))
	overlap_ratio <-  round(bin_width / bin_change)

	NN <- nbins * N.cells

	first_spike <- min(which(spikes[,"time"] >= bin.start.end[1,1]))
	last_spike <- max(which(spikes[,"time"] < bin.start.end[nbins,2]))
	spt_pos <- spikes[first_spike:last_spike,, drop=F]

	N.spikes <- nrow(spt_pos)
	
	# sparseness <- N.spikes * sum(bin.start.end[,2] - bin.start.end[,1]) / (max(bin.start.end) - min(bin.start.end)) / NN
	# if (sparseness < 1/50)	rast <- Matrix(0, nbins, N.cells)	else rast <- matrix(0, nbins, N.cells)
	rast <- matrix(0, nbins, N.cells)
		
	## this is a fairly quick algorithm - but only works for ordered and continuous bins
	i_start <- 1	
	i_end <- min(nbins, i_start + 10)
	for (i in 1:N.spikes){
		t.sp <- spt_pos[i,"time"]
		
		bin_not_found <- T
		while(bin_not_found){
			i_ends_later <- which(bin.start.end[i_start:i_end,2] > t.sp)
			j.starts_earlier <- which(bin.start.end[i_start:i_end,1] <= t.sp)
			if (length(i_ends_later) > 0){
				i.bin <- i_start -1 + min(i_ends_later) # first bin: first that ends later				
				j.bin <- i_start -1 + max(j.starts_earlier) # last bin: last that starts earlier
				if ((min(i_ends_later) > 0) & (max(j.starts_earlier) < (i_end - i_start + 1))) {
					bin_not_found <- F 
				} else if ((min(i_ends_later) > 0) & (max(j.starts_earlier) == (i_end - i_start + 1)) & (i_end==nbins)) {
					bin_not_found <- F 
				} else {
					i_end <- min(nbins, i_end + 10)
				}
			}  else {
				i_end <- min(nbins, i_end + 10)
			}
		}
		# if (t.sp < bin.start.end[i.bin,2]){# we don't need to check this if the bins are increasing...
		if (renumber_cells) i.cell <- which( ii.cells == spt_pos[i,"cell"]) else i.cell <- spt_pos[i, "cell"]
		rast[i.bin:j.bin, i.cell] <- rast[i.bin:j.bin, i.cell] + 1
		# }
		i_start <- i.bin
		i_end <- min(nbins, i_start + 10)
	}
	
	if (sum(rast) > overlap_ratio*nrow(spt_pos)) stop('more spikes in raster than expected...')
	if (sum(rast) < (nbins - 4 * overlap_ratio) / nbins * overlap_ratio*nrow(spt_pos)) stop('too few spikes in raster compared to input spike train ...')

	# image(rowMeans(bin.start.end[,1:2]), 1:N.cells, rast)
	# points(spikes[,1], spikes[,2], pch=1, col=grey(0.5, alpha=0.25))
	rast
}

sp2thetaRasts <- function(spikes, theta_chunks, N_bins=3, N.cells){
	# a function to separate spikes of a train according to their theta phase.
	# return value: a list with multiple rasters, each corresponding to a particular theta phase
	# each theta cycle is divided to three bins of equal number of spikes. 
	# the bin width is also returned
	
	if (N_bins != 3) stop('N_bins should be 3!')
	
	n_periods <- nrow(theta_chunks)
	ispikes_in_theta <- ((spikes[,1] >= theta_chunks[1,1]) & (spikes[,1] < theta_chunks[n_periods,2]))
	spikes <- spikes[ispikes_in_theta, ]

	N_sp <- nrow(spikes)
	sp_ThetaCycle <- rep(NA, N_sp)
	sp_ThetaPhase <- rep(NA, N_sp)
	chunk_mids <- matrix(NA, n_periods, N_bins-1)
	for (i_theta in 1:n_periods){
		ispikes_in_bin <- which((spikes[,1] >= theta_chunks[i_theta,1]) & (spikes[,1] < theta_chunks[i_theta,2]))
		spikes_i <- matrix(spikes[ispikes_in_bin,], ncol=2)
		n <- length(ispikes_in_bin)
		if (n >= N_bins){
			sp_phase <- sort(c(rep(1:N_bins, floor(n / N_bins)), sample(1:N_bins, N_bins))[1:n])
			sp_ThetaPhase[ispikes_in_bin] <- sp_phase
			sp_ThetaCycle[ispikes_in_bin] <- i_theta
			for (i_phase in 2:N_bins){
				t_start_new <- spikes_i[min(which(sp_phase == i_phase)),1]
				t_end_old <- spikes_i[max(which(sp_phase == i_phase-1)),1]
				chunk_mids[i_theta, i_phase-1] <- (t_end_old + t_start_new) / 2
			}
		} else { # 
			sp_phase <- rep(1, n)
			for (i_phase in 2:N_bins){
				chunk_mids[i_theta, i_phase-1] <- theta_chunks[i_theta,1]  +  (i_phase - 1) * (theta_chunks[i_theta,2] - theta_chunks[i_theta,1]) / N_bins
				sp_phase[spikes_i[,1] > chunk_mids[i_theta, i_phase-1]] <- i_phase
			}
			sp_ThetaPhase[ispikes_in_bin] <- sp_phase
			sp_ThetaCycle[ispikes_in_bin] <- i_theta			
		}
	}
	
	# valid_spikes <- spikes[!is.na(sp_ThetaCycle),]
	# valid_spPhase <- sp_ThetaPhase[!is.na(sp_ThetaCycle)]
	spikes <- spikes[!is.na(sp_ThetaPhase),]
	sp_ThetaPhase <- sp_ThetaPhase[!is.na(sp_ThetaPhase)]
	
	rast1 <- spt2raster_bins(spikes[sp_ThetaPhase ==1,,drop=F], theta_chunks, N.cells=N.cells)
	bin_size1 <- chunk_mids[,1] - theta_chunks[,1]
	bin_mids1 <- (chunk_mids[,1] + theta_chunks[,1])/2
	bintime1 <- cbind(bin_size1, bin_mids1)
	colnames(bintime1) <- c('duration', 'time')
	
	rast2 <- spt2raster_bins(spikes[sp_ThetaPhase ==2,,drop=F], theta_chunks, N.cells=N.cells)
	bin_size2 <- chunk_mids[,2] - chunk_mids[,1]
	bin_mids2 <- (chunk_mids[,2] + chunk_mids[,1])/2
	bintime2 <- cbind(bin_size2, bin_mids2)
	colnames(bintime2) <- c('duration', 'time')

	rast3 <- spt2raster_bins(spikes[sp_ThetaPhase==3,,drop=F], theta_chunks, N.cells=N.cells)
	bin_size3 <- theta_chunks[,2] - chunk_mids[,2]
	bin_mids3 <- (theta_chunks[,2] + chunk_mids[,2])/2
	bintime3 <- cbind(bin_size3, bin_mids3)
	colnames(bintime3) <- c('duration', 'time')
	
	theta_rasters <- list(raster1=rast1, raster2=rast2, raster3=rast3, bintime1=bintime1, bintime2=bintime2, bintime3=bintime3)
	theta_rasters
}

################################################################################
## posterior mean decoding using ratemaps
## this function calculates to posterior mean position in each timestep using the spikes of the neurons and the prior over the positions
### inputs: spikes - vector with the number of spikes for each of the N cells 
## 	ratemaps: an array of N x K x K dimansions with the firing rates 
## 	posmap: a K x K array with the time spent at each location on a K x K grid. Its attributes 'xcenters' and 'ycenters' provide the x and y coordinates
##		
## output: a list with the following elements: 
## 	estim_pos: posterior mean position
##		estim_var: posterior variance

postmean_decode <- function(spikes, ratemaps, posmap, dt=0.1, graphics=F){

	prior <- as.vector(posmap / sum(posmap))
	xx <- attr(posmap, 'xcenters')
	yy <- attr(posmap, 'ycenters')
	nx <- length(xx)
	ny <- length(yy)
	xmat <- matrix(rep(xx, each=ny), nx, byrow=T)
	ymat <- matrix(rep(yy, nx), nx, byrow=T)
		
	min_rate <- 1 / sum(posmap)
	ratemaps[ratemaps < min_rate] <- min_rate / 2

	lambda_dt <- ratemaps * dt
	n.cells <- dim(ratemaps)[1]
	n.pixels <- prod(dim(ratemaps)[2:3])
	dim(lambda_dt) <- c(n.cells, n.pixels)	

	NAcells <- is.na(spikes)
	valid_cells <- !NAcells
	
	# rerate <- lambda_dt
	# dim(rerate) <- c(n.cells, dim(ratemaps)[2], dim(ratemaps)[3])
	prior2 <- log(prior) - colSums(lambda_dt[valid_cells,])
	
	LL <- prior2 + as.vector(spikes[valid_cells] %*% log(lambda_dt[valid_cells,]))
	LL <- array(LL, dim=c(dim(ratemaps)[2], dim(ratemaps)[3]))
	LL <- normaliseLL(LL)

	post <- exp(LL)
	post <- post / sum(post)
	
	xhat <- sum(xmat * post)
	yhat <- sum(ymat * post)
	# xhat <- xmat[which.max(post)]
	# yhat <- ymat[which.max(post)]
	estim_pos <- c(xhat, yhat)
	estim_var <- c(sum(xmat^2 * post) - xhat^2, sum(ymat^2 * post) - yhat^2)
	if (graphics > 0){
		image(xx, yy, exp(LL), col=colormap(colormaps$hot, 48))
		points(xhat, yhat, pch=21, col=rgb(0, 1, 1), bg=rgb(0, .5, .5))
	}

	decode <- list(LL=LL, mean=estim_pos, var=estim_var)
	decode
}

################################################################################
## posterior mean decoding using ratemaps
## this function calculates to posterior mean position in each timestep using the spikes of the neurons and the prior over the positions
### inputs: chunks: list of running epochs. Each element contains a list with the spike count matrix ('rast', a T [time steps] x N [# of cells] matrix) and the true position ('pos', T x 2 matrix). 
## 	ratemaps: an array of N x K x K dimansions with the firing rates 
## 	posmap: a K x K array with the time spent at each location on a K x K grid. Its attributes 'xcenters' and 'ycenters' provide the x and y coordinates
##		
## output: a list with lists corresponding to the chunks
## each with the following elements: 
## 	estim_pos: posterior mean position
##		estim_var: posterior variance
## 	rsqerr: root mean squared error
## 	nspikes: number of spikes

# decode_ratemaps <- function(chunks, ratemaps, posmap, dt=0.1, graphics=F, nshift=0){
postmean_chunks <- function(chunks, ratemaps, posmap, dt=0.1, graphics=F, nshift=0){
	n.chunks <- length(chunks)
	
	i.chunk <- 1
	prior <- as.vector(posmap / sum(posmap))
	xx <- attr(posmap, 'xcenters')
	yy <- attr(posmap, 'ycenters')
	nx <- length(xx)
	ny <- length(yy)
	xmat <- matrix(rep(xx, each=ny), nx, byrow=T)
	ymat <- matrix(rep(yy, nx), nx, byrow=T)
		
	min_rate <- 1 / sum(posmap)
	ratemaps[ratemaps < min_rate] <- min_rate / 2

	lambda_dt <- ratemaps * dt
	n.cells <- dim(ratemaps)[1]
	n.pixels <- prod(dim(ratemaps)[2:3])
	dim(lambda_dt) <- c(n.cells, n.pixels)	

	NAcells <- is.na(chunks[[1]]$rast[1,])
	valid_cells <- !NAcells
	
	# rerate <- lambda_dt
	# dim(rerate) <- c(n.cells, dim(ratemaps)[2], dim(ratemaps)[3])
	prior2 <- log(prior) - colSums(lambda_dt[valid_cells,])
	decode <- list()	
	
	for (i.chunk in 1:n.chunks){
		section <- chunks[[i.chunk]]
		n.steps <- nrow(section$rast)
		
		estim_pos <- matrix(NA, n.steps, 2)
		estim_var <- matrix(NA, n.steps, 2)
		rsqerror <- matrix(NA, n.steps, 2*nshift+1)
		
		for (i.step in 1:n.steps){
			code <- section$rast[i.step,]
			LL <- prior2 + as.vector(code[valid_cells] %*% log(lambda_dt[valid_cells,]))
			LL <- array(LL, dim=c(dim(ratemaps)[2], dim(ratemaps)[3]))
			LL <- normaliseLL(LL)

			if (graphics > 0){
				image(xx, yy, exp(LL), col=colormap(colormaps$hot, 48))
				points(section$pos[i.step,1], section$pos[i.step,2], pch=21, col=1, bg=7)
			}
						

			post <- exp(LL)
			post <- post / sum(post)
			
			xhat <- sum(xmat * post)
			yhat <- sum(ymat * post)
			# xhat <- xmat[which.max(post)]
			# yhat <- ymat[which.max(post)]
			estim_pos[i.step,] <- c(xhat, yhat)
			estim_var[i.step,] <- c(sum(xmat^2 * post) - xhat^2, sum(ymat^2 * post) - yhat^2)
			for (k in (1:(2*nshift+1))){
				kk <- k - nshift - 1
				if (((i.step+kk) > 0) & ((i.step+kk) <= n.steps)){
					rsqerror[i.step, k] <- sqrt(sum((c(xhat, yhat) - section$pos[i.step+kk,])^2))	
				}
			}	
		}
		
		# plot(section$pos[,1], section$pos[,2], col=viridis(n.steps), pch=16, xlim=c(0, 200), ylim=c(0, 200))
		# points(estim_pos[,1], estim_pos[,2], col=viridis(n.steps), pch=17, xlim=c(0, 200), ylim=c(0, 200))

		decode[[i.chunk]] <- list(estim_pos=estim_pos, estim_var=estim_var, rsqerror=rsqerror, nspikes=rowSums(section$rast, na.rm=T))
		cat(i.chunk, ' ')			
	}
	decode
}

####################################################
## filtering of the observed spikes - assumes Gaussian random walk dynamics with 
dynamic_decode <- function(log_prior, rast, ratemaps, xx, sigma=3, dt=0.005){
	# log_prior: log posterior of the previous timestep
	# rast: raster of L x N with the spike counts for each cell in each timestep
	# ratemaps: an array of N x K x K dimansions with the firing rates 
	# xx: x and y coordinates of the ratemaps and the log_prior - they should be square matrices
	# sigma: spatial variance of the random walk in sigma^2 is in cm^2/s
	# dt: temporal resolution (s) can not be very small, since sigma needs to be at least as large as the resolution of the ratemaps 
	# the SD of the Gaussian blur will be ss = sigma / xw * sqrt(deltaT)
	#
	# P(x_t | s_1:t) = P(x_t | x_t-1, s_1:t-1) P(s_t|x_t)
	
	L <- nrow(rast)
	K <- dim(ratemaps)[2]

	lambda_dt <- ratemaps * dt
	n.cells <- dim(ratemaps)[1]
	n.pixels <- K*K
	dim(lambda_dt) <- c(n.cells, n.pixels)	

	NAcells <- is.na(rast[1,])
	valid_cells <- !NAcells

	nx <- length(xx)
	xmat <- matrix(rep(xx, each=nx), nx, byrow=T)
	ymat <- matrix(rep(xx, nx), nx, byrow=T)

	xw <- (xx[2] - xx[1]) * 100 # cm
	deltaT <- 0.005 # ms
	ss <- sigma / xw * sqrt(deltaT*1000) # i


	if (L < 100) post <- array(NA, dim=c(L, K, K))
	means <- matrix(NA, L, 2)
	vars <- matrix(NA, L, 2)
	
	for (ell in 1:L){
		prior <- exp(log_prior)
		pred <- 	blur(as.im(prior), ss)$v
		i_negative <- which(pred < 0)
		pred[i_negative] <- prior[i_negative]
		
		Lpost <- log(pred) - colSums(lambda_dt[valid_cells,]) + as.vector(rast[ell,valid_cells] %*% log(lambda_dt[valid_cells,]))
		Lpost <- array(Lpost, dim=c(K, K))
		Lpost <- normaliseLL(Lpost)
		log_prior <- Lpost
		
		epost <- exp(Lpost)
		epost <- epost / sum(epost)

		xhat <- sum(xmat * epost)
		yhat <- sum(ymat * epost)
		means[ell,] <- c(xhat, yhat)
		vars[ell,] <- c(sum(xmat^2 * epost) - xhat^2, sum(ymat^2 * epost) - yhat^2)
		if (L < 100) post[,,ell] <- epost
	}
	
	out <- list(post=post, means=means, vars=vars)
	if (L < 100) out$post <- post
	out
}


# map1 <- matrix(0, 40, 40)
# map1[20, 20] <- 1

# xw <- 5 # cm
# sigma <- 3 # cm
# deltaT <- 1 # ms
# ss <- sigma / xw * sqrt(deltaT)
# map2 <- blur(as.im(map1), ss)$v
# for ( i in 1:4) map2 <- blur(as.im(map2), ss)$v

# deltaT <- 5 # ms
# ss <- sigma / xw * sqrt(deltaT)
# map3 <- blur(as.im(map1), ss)$v

# range(map2-map3)
# image(map3)
# image(map2)


################################################################################
## normalising log likelihoods
## it is often a problem, that log likelioods get very small or very large, so when exponentiating them, we get NaN
## in practice, LL 
## - should be between -700 and 700
## - LL is invariant for additions (unnormalised)
## - small values are negligible

normaliseLL <- function(LL){
	if (max(LL) < -500) LL <- LL - max(LL) # we make it larger
	if (max(LL) > 500) LL <- LL - max(LL) + 100 # we make it smaller
	LL
}


################################################################################
## DDC likelihood 1D - in a form where it is relatively easy to optimise the parameters for ML decoding

DDC_likelihood_1D <- function(pars, sp, rates_vec, x.vec, deltaT, opars=NULL){
	## sampling the encoded distribution on a regular grid
	if ('mu' %in% names(pars)) mu <- pars['mu'] else mu <- opars$mu
	if ('sigma' %in% names(pars)) sigma <- pars['sigma'] else sigma <- opars$sigma
	
	P.encoded <- dnorm(x.vec, mu, sigma) # prepare the distribution: a 1D function
	P.encoded <- P.encoded / sum(P.encoded)
	
	## encoding rates via convolution of the distribution with the receptive field
	rate <- P.encoded %*% rates_vec * deltaT
	cellids <- which(!is.na(sp))
	LL <- sum(dpois(sp[cellids], rate[cellids], log=T))
  	LL
}

 ################################################################################
## DDC likelihood - in a form where it is relatively easy to optimise the parameters for ML decoding

DDC_likelihood <- function(pars, sp, rates_vec, g.x.vec, g.y.vec, deltaT, opars=NULL){
  
	## sampling the encoded distribution on a regular grid
	if ('mu1' %in% names(pars)) mu1 <- pars['mu1'] else mu1 <- opars$mu1
	if ('mu2' %in% names(pars)) mu2 <- pars['mu2'] else mu2 <- opars$mu2
	if ('sigma' %in% names(pars)) sigma <- pars['sigma'] else sigma <- opars$sigma
	
	pars.encoded <- list(mu=matrix(c(mu1, mu2), 1), sds=matrix(c(sigma, sigma), 1), rmax=1, r0=0) # extract the parameters to be encoded
	P.encoded <- get.rate.x(rbind(g.x.vec, g.y.vec), pars.encoded)	# prepare the distribution: a 2D function
	P.encoded <- P.encoded / sum(P.encoded)
	
	## encoding rates via convolution of the distribution with the receptive field
	rate <- P.encoded %*% rates_vec * deltaT
	cellids <- which(!is.na(sp))
	LL <- sum(dpois(sp[cellids], rate[cellids], log=T))
  	LL
}

################################################################################
## calculates the posterior mean and variance of the position given the spikes
## this function calculates to posterior mean position in each timestep using the spikes of the neurons and the prior over the positions
### inputs: chunks: list of running epochs. Each element contains a list with the spike count matrix ('rast', a T [time steps] x N [# of cells] matrix) and the time, true position, head direction and bin width ('pos', T x 5 matrix). 
## 	ratemaps: an array of N x K x K dimansions with the firing rates 
## 	posmap: a K x K array with the time spent at each location on a K x K grid. Its attributes 'xcenters' and 'ycenters' provide the x and y coordinates
##		
## output: a list with lists corresponding to the chunks
## each with the following elements: 
## 	estim_pos: posterior mean position
##		estim_var: posterior variance
## 	rsqerr: root mean squared error
## 	nspikes: number of spikes

decode_ratemaps <- function(raster, bin_widths, ratemaps, posmap, MAP=FALSE, dim=2){
	prior <- as.vector(posmap / sum(posmap))
	min_rate <- 1 / sum(posmap)
	ratemaps[ratemaps < min_rate] <- min_rate / 2
	n.cells <- dim(ratemaps)[1]
	NAcells <- is.na(raster[1,])
	valid_cells <- !NAcells

	xx <- attr(posmap, 'xcenters')
	nx <- length(xx)
	if (dim == 2){
		yy <- attr(posmap, 'ycenters')
		ny <- length(yy)
		xmat <- matrix(rep(xx, each=ny), nx, byrow=T)
		ymat <- matrix(rep(yy, nx), nx, byrow=T)
		n.pixels <- prod(dim(ratemaps)[2:3])
	} else {
		n.pixels <- dim(ratemaps)[2]		
	}
		
	n.steps <- nrow(raster)
	if (!is.null(n.steps)){
		estim_pos <- matrix(NA, n.steps, dim)
		estim_var <- matrix(NA, n.steps, dim)
		
		for (i.step in 1:n.steps){
			code <- raster[i.step,]
	
			lambda_dt <- ratemaps * bin_widths[i.step]
			if (dim == 2){
				dim(lambda_dt) <- c(n.cells, n.pixels)
				prior2 <- log(prior) - colSums(lambda_dt[valid_cells,])	
				LL <- prior2 + as.vector(code[valid_cells] %*% log(lambda_dt[valid_cells,]))
				LL <- array(LL, dim=c(dim(ratemaps)[2], dim(ratemaps)[3]))
				LL <- normaliseLL(LL)
				post <- exp(LL)
				post <- post / sum(post)
						
	
				if (MAP){
					xhat <- xmat[which.max(post)]
					yhat <- ymat[which.max(post)]
				} else {
					xhat <- sum(xmat * post) # posterior mean X
					yhat <- sum(ymat * post) # posterior mean Y				
				}
				estim_pos[i.step,] <- c(xhat, yhat)
				estim_var[i.step,] <- c(sum(xmat^2 * post) - xhat^2, sum(ymat^2 * post) - yhat^2)
			} else {
				prior2 <- log(prior) - colSums(lambda_dt[valid_cells,])	
				LL <- prior2 + as.vector(code[valid_cells] %*% log(lambda_dt[valid_cells,]))
				LL <- normaliseLL(LL)
				post <- exp(LL)
				post <- post / sum(post)
							
				if (MAP){
					xhat <- xx[which.max(post)]
				} else {
					xhat <- sum(xx * post) # posterior mean X
				}
				estim_pos[i.step,] <- xhat
				estim_var[i.step,] <- sum(xx^2 * post) - xhat^2		
			}				
		}
	}
	post <- list(mean=estim_pos, var=estim_var)
	post
}
################################################################################
## similar to decode ratemaps, but 
## - runs on theta chunks - periods with ongoing theta separated by some time interval when the animal stopped
## - also calculates the error wrt the temporally or spatially shifted x and y coordinates
## - 
## this function calculates to posterior mean position in each timestep using the spikes of the neurons and the prior over the positions
### inputs: chunks: list of running epochs. Each element contains a list with the spike count matrix ('rast', a T [time steps] x N [# of cells] matrix) and the time, true position, head direction and bin width ('pos', T x 5 matrix). 
## 	ratemaps: an array of N x K x K dimansions with the firing rates 
## 	posmap: a K x K array with the time spent at each location on a K x K grid. Its attributes 'xcenters' and 'ycenters' provide the x and y coordinates
##		
## output: a list with lists corresponding to the chunks
## each with the following elements: 
## 	estim_pos: posterior mean position
##		estim_var: posterior variance
## 	rsqerr: root mean squared error
## 	estim_dir: deviation between the estimated and the actual position relative to the heading of the animal

decode_theta_chunks <- function(chunks, ratemaps, posmap, graphics=F, calc_error=F, tshift=0, xshift=0, txyv=NULL, MAP=FALSE, dim=2){

	if (length(tshift) != length(xshift)) {
		print('tshift and xshift have different length! We set both to 0.')
		tshift <- 0 # time shift (s) used to calculate the decoding error. The reference position taken at t+tshift.
		xshift <- 0 # statial shift (cm) used to calculate the decoding error. The reference position taken as x + xshift in the running direction.
	}
	if (sum((tshift != 0) & (xshift!=0))){
		print('tshift and xshift have elements where both of them are non-zero! We set both to 0.')
		tshift <- 0 # time shift (s) used to calculate the decoding error. The reference position taken at t+tshift.
		xshift <- 0 # statial shift (cm) used to calculate the decoding error. The reference position taken as x + xshift in the running 	
	}
	
	if (max(abs(tshift)) != 0){
		if (is.null(txyv)){
			stop ('txyv must be provided if tshift is not 0!')
		}
	}
	n.shifts <- length(tshift)
	n.chunks <- length(chunks)
		
	decode <- list()	
	
	xx <- attr(posmap, 'xcenters')
	yy <- attr(posmap, 'ycenters')
	
	for (i.chunk in 1:n.chunks){
		section <- chunks[[i.chunk]]
		n.steps <- nrow(section$rast)
		if (!is.null(n.steps)){
			post <- decode_ratemaps(section$rast, bin_widths=section$pos[,'bin_width'], ratemaps, posmap, MAP=MAP, dim=dim)
			decode[[i.chunk]] <- list(estim_pos=post$mean, estim_var=post$var)
	
			if (calc_error){
				estim_pos <- post$mean
				rsqerror <- matrix(NA, n.steps, n.shifts)
				estim_dir <- matrix(NA, n.steps, n.shifts)
					
				for (i.step in 1:n.steps){
					if (dim==2){
						xhat <- estim_pos[i.step,1]
						yhat <- estim_pos[i.step,2]
						motion_dir <- section$pos[i.step,4]
		
						for (i.shift in 1:n.shifts){
							xsh <- xshift[i.shift]
							tsh <- tshift[i.shift]					
							# ref_pos <- section$pos[i.step,2:3] + c(xshift*cos(motion_dir), xshift*sin(motion_dir))
							if (tsh == 0){
								ref_pos <- section$pos[i.step,2:3] + c(xsh*cos(motion_dir), xsh*sin(motion_dir))	
							} else {
								ref_time <- section$pos[i.step,1] + tsh
								ref_pos <- c(approx(txyv[,1], txyv[,2], ref_time)$y, approx(txyv[,1], txyv[,3], ref_time)$y)
							}
						
							abs_estim_dir <- atan2(yhat - ref_pos[2], xhat - ref_pos[1])
							estim_dir[i.step, i.shift] <- (abs_estim_dir - motion_dir + pi) %% (2*pi) - pi										
							rsqerror[i.step, i.shift] <- sqrt(sum((c(xhat, yhat) - ref_pos)^2))
						}
		
						if (graphics > 0){
							image(xx, yy, posmap, col=colormap(colormaps$hot, 48))
							points(section$pos[i.step,2], section$pos[i.step,3], pch=21, col=1, bg=7)
							points(xhat, yhat, pch=21, col=1, bg=viridis(10)[7])
						} 
					} else {
						motion_dir <- sign(section$pos[n.steps,2] - section$pos[1,2])
						xhat <- estim_pos[i.step,1]
						for (i.shift in 1:n.shifts){
							xsh <- xshift[i.shift]
							tsh <- tshift[i.shift]					
							if (tsh == 0){
								ref_pos <- section$pos[i.step,2] + xsh*motion_dir
							} else {
								ref_time <- section$pos[i.step,1] + tsh
								ref_pos <- approx(txyv[,1], txyv[,2], ref_time)$y
							}
							rsqerror[i.step, i.shift] <- xhat - ref_pos
						}
					}
				}
				decode[[i.chunk]]$rsqerror <- rsqerror
				if (dim==2) decode[[i.chunk]]$dir <- estim_dir
			}
		} else {
			decode[[i.chunk]] <- list(estim_pos=NA, estim_var=NA, rsqerror=NA, dir=NA)			
		}
	}
	decode
}



################################################################################
## fit_Gauss_ratemap(spi, pos, vv, posmap.raw, dt)

# pos <- txyv


fit_Gauss_ratemap <- function(spi, txy, i_data, posmap, graphics=0){
	## i_data: index for the position to include in the estimate
	# get the spikes with velocity over threshold

	dt <- median(diff(txyv[,1]))

	ii_run <- approx(txy[,1], i_data, spi, method='constant')$y
	spi_run <- spi[!!ii_run]

	x_spikes <- approx(txy[,1], txy[,2], spi_run)$y
	y_spikes <- approx(txy[,1], txy[,3], spi_run)$y
	xyn <- cbind(x_spikes, y_spikes, rep(1, length(x_spikes)))
		
	pars <-c(mu1=100, mu2=100, sd1=50, sd2=50, alpha=-1, beta=0.1*dt) # params are in 1/(33 s) - so not in Hz!
	rate <- nrow(xyn) / (sum(i_data) * dt) # in Hz
	pars['alpha'] <- log(rate*dt * 40000 / (2*pi*sqrt(det(diag(c(pars['sd1']^2, pars['sd2']^2))))))

	if (rate > 5) {
		print('approximating LL based on binned spike counts')
		ratemap.raw <- estimate_ratemap(spi, txy, i_data, xw=10, yw=10, sigma=1, posmap=posmap, graphics=0)
		counts <- ratemap.raw * posmap
		opars <- optim(pars, Gauss_LL_map, gr=grad_Gauss_LL_map, xyn=xyn, posmap=posmap, dt=dt, countmap=counts, method='BFGS', control=list(trace=5, maxit=100, fnscale=-1))
	} else {
		opars <- optim(pars, Gauss_LL_map, gr=grad_Gauss_LL_map, xyn=xyn, posmap=posmap, dt=dt, method='BFGS', control=list(trace=5, maxit=100, fnscale=-1))	
	}


	if (graphics > 0){	
		rr_map <- matrix(NA, 20, 20)
		for (i in 1:20){
			for (j in 1:20){
				x <- i*10-5
				y <- j*10-5
				rr_map[i,j] <- Gauss_rate(opars$par, c(x,y))
				# rr_map[i,j] <- Gauss_rate(pars, c(x,y))
			}
		}

		maxr <- ceiling(max(rr_map)/dt * 1.1)
		minr <- 0
		
		image(1:20*10-5, 1:20*10-5, rr_map/dt, col=colormap(colormaps$hot, 48), xlab='', ylab='', axes=F, br=seq(minr, maxr, length=49))		
		if (rr_map[3, 18]/dt / maxr < 0.5) textcol=grey(1) else textcol=grey(0)
		if (graphics > 1)	points(x_spikes, y_spikes, pch=16, cex=0.4, col=viridis(6, option='D', alpha=0.5)[4])
		text(25, 175, paste(maxr, 'Hz'), col=textcol)		
	}
	opars
}

##########################################################################
ratemap_LL <- function(posmap, ratemap, countmap, dt){
	# LL (spikes | ratemap)
	LL <- 0
	for (i in 1:nrow(posmap)){
		for (j in 1:ncol(posmap)){
			lambda <- ratemap[i, j]
			# LL <- LL - lambda * dt * posmap[i, j] / dt + log(lambda) * countmap[i, j]
			LL <- LL - lambda * posmap[i, j] + log(lambda*dt) * countmap[i, j]
		}
	}
	LL
}


##########################################################################
skaggs93.info <- function(r.x, Px){
	## I = \sum_x r(x)p(x) log(r(x)/r)
	## r.x: estimated firing rate as function of x (spikes/sec)
	## Px: P(x), unnormalised
	## output: information rate, bit/s
	## Px <- total.t; r.x <- rate.x[,43]
	Px <- Px / sum(Px)
	r.x[r.x < 1e-10] <- 1e-10
	r <- sum(Px * r.x)
	KL <- r.x * Px * log(r.x/r, 2)
	KL[r.x==0] <- 0
	I <- sum(KL)
	I
}

##########################################################################
Fisher_error <- function(rates, dx, Tbin){
	## calculate the lower bound on estimation error based on Fisher info in 2D place cells
	## rates: ratemaps of the N cells
	## dx: the x resolution (e.g., in cm)
	## Tbin the duration of the observation interval (in seconds)
	##
	## output: 
	## lower bound on the estimation error (in the same units as dx, e.g, in cm)
	
	N <- dim(rates)[1]
	K <- dim(rates)[2]
	L <- dim(rates)[3]
		
	rates2 <- array(NA, dim=c(N, K-1, L-1))
	for (k in 1:(K-1)){
		for (l in 1:(L-1)){
			rates2[,k, l] <- (rates[,k, l] + rates[,k+1, l] + rates[,k, l+1] + rates[,k+1, l+1]) / 4
		}
	}
	
	slopes <- array(NA, dim=c(2, N, K-1, L-1))
	for (n in 1:N){
		slope <- get_slope(rates[n,,], dx=dx)
		slopes[1, n,,] <- slope$dx
		slopes[2, n,,] <- slope$dy		
	}
	
	IFx <- Tbin * apply(slopes[1,,,]^2/rates2, c(2,3), sum)
	IFy <- Tbin * apply(slopes[2,,,]^2/rates2, c(2,3), sum)
	errFisher <- sqrt(1/IFx + 1/IFy) ## cm
	# cm = sqrt(1/dT*slope)
	# cm^2 = 1/dT * slope
	errFisher
}



## A function to numerically estimate the slope of an image
## If the image has dimensions (N+1, N+1) then the output is 2 (N, N) matrices, the partial derivaties in the x and y direction
## so the estimation is done for the vertices of hte matrix
##
## We used this function to estimate the Fisher information of a population

get_slope <- function(mat, dx){
# input: mat: N x M matrix, representing the value of a function F(x) on regular grids with 
# 			dx: the distance between grid points
# output: list of two (N-1) x (M-1) matrices with the numerical estimate of the partial derivatives dF/dx and dF/dy
# 			if mat was given at points x1 and x2 then the output is calculated for (x1 + x2) / 2
	N <- dim(mat)[1]
	M <- dim(mat)[2]
	dmat.dx <- matrix(NA, N-1, M)
	dmat.dy <- matrix(NA, N, M-1)
	for (x in 1:(N-1)) dmat.dx[x,] <- (mat[x+1,] - mat[x,]) / dx
	for (y in 1:(M-1)) dmat.dy[,y] <- (mat[,y+1] - mat[,y]) / dx

	dmat.dxx <- matrix(NA, N-1, M-1)
	dmat.dyy <- matrix(NA, N-1, M-1)
	for (x in 1:(N-1)) dmat.dyy[x,] <- (dmat.dy[x,] + dmat.dy[x+1,]) / 2
	for (y in 1:(M-1)) dmat.dxx[,y] <- (dmat.dx[,y] + dmat.dx[,y+1]) / 2
	
	slope <- list(dx = dmat.dxx, dy = dmat.dyy)
	slope	
}



####################################################
## apply shuffling in theta chunks
## INPUT:
## 	splist: a list with the theta chunks
## 	each theta chunk has a raster (matrix of N_bins x N_cells) with the spike counts
## OUTPUT:
## 	splist with shuffled spikes

shuffle_spikes <- function(splist){
	allrast <- splist[[1]]$rast
	for (j in 2:length(splist)){
		allrast <- rbind(allrast, splist[[j]]$rast)
	}
	Ncells <- ncol(allrast)
	allrast <- allrast[which(rowSums(is.na(allrast))<Ncells),]
	
	shrast <- matrix(0, nrow(allrast), ncol(allrast))
	for (i_cell in 1:ncol(allrast)){
		shrast[,i_cell] <- permute(allrast[,i_cell])
	}

	k1 <- 1
	k2 <- 0
	for (j in 1:length(splist)){
		if (!is.na(splist[[j]]$pos[1,1])){
			k2 <- k2 + nrow(splist[[j]]$rast)
			splist[[j]]$rast <- shrast[k1:k2,]
			k1 <- k2 + 1
		} else {
			splist[[j]]$rast <- NA
		}
	}
	splist
}


########################################################
calc_map_SD <- function(ratemaps, x){
	# ratemaps: each row is the firing rate of a cell
	N <- nrow(ratemaps)
	SDs <- rep(NA, N)
	Ms <- rep(NA, N)
	for (i in 1:N) {
		PP <- ratemaps[i,] / sum(ratemaps[i,])	
		mm <- PP %*% x
		ss <- PP %*% x^2
		SDs[i] <- ss - mm^2 
		Ms[i] <- mm
	}
	cbind(Ms, SDs)
}

