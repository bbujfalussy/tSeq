min.err.decoded.true <- function(xy.decoded, xy.true){
	## calculates the error between the decoded trajectory and the true trajectory
	## by finding the closest apposition between the true trajectory and the decoded points
	## x and y coordinates are the rows of the matrices
	nseg <- ncol(xy.decoded)
	err <- rep(NA, nseg)
	for (t in 1:nseg){
		xy.t <- xy.decoded[,t]
		dist.points <- sqrt(colSums(xy.true - xy.t)^2)
		err[t] <- min(dist.points)
	}
	err
}

err.closest <- function(xy.decoded, xy.true){
	## xy.true: the true trajectory - for the entire experiment, a matrix of 3 x N, with the rows t, x, y
	## xy.decoded: the decoded trajectory - for each theta cycle, an array of 2 x Nseg X N
	## calculates the error between the decoded trajectory and the true trajectory
	## by finding the closest apposition between the true trajectory and the decoded points
	N <- ncol(xy.true)
	Nseg <- dim(xy.decoded)[2]
	err <- matrix(NA, Nseg, N)
	for (t in 11:(N-10)){
		xy.t <- xy.true[2:3,(t-10):(t+10)]
		xy.d.t <- xy.decoded[,,t]
		err[,t] <- min.err.decoded.true(xy.d.t, xy.t)		
	}
	err
}

### a function to prepare a spike count vector
spt2spc <- function(spt, N.cells){
	spc <- rep(0, N.cells)
	for (i in 1:length(spt)) spc[spt[i]] <- spc[spt[i]] + 1
	spc
}


## this function plots the mean of a set of samples and a shaded interval around them at mean ± nsd * sd.
polyplot.sam <- function(t=NULL, samples, col.bg=grey(0.8), col=1, xlab=NULL, ylab=NULL, ylim=NULL, xlim=NULL, add=F, axes=T, tit=NULL, nsd=1){


	if (is.null(t)) t <- seq(1, ncol(samples))
	
	m.mat <- colMeans(samples, na.rm=T)
	sd.mat <- nsd * apply(samples, 2, sd, na.rm=T)

	tt <- c(t, rev(t))
	ss <- c(m.mat+sd.mat, rev(m.mat-sd.mat))
	if (is.null(ylim)) ylim <- range(ss)
	if (is.null(xlim)) xlim <- range(t)
	if (add==F)	plot(t, m.mat, ylim=ylim, xlim=xlim, t="n", axes=F, xlab=xlab, ylab=ylab, main=tit)
	polygon(tt,ss, col=col.bg, border=grey(1))
	if ((!add) & (axes)) {axis(1); axis(2, las=2)}
	lines(t, m.mat, lwd=2, col=col)
	lines(t, m.mat+sd.mat, lwd=1, lty=2, col=col)
	lines(t, m.mat-sd.mat, lwd=1, lty=2, col=col)

}

##########################################
## adding scalebar to plots


scalebar <- function(x1, x2, y1, y2, xscale=NULL, yscale=NULL){
	lines(c(x1, x1, x2), c(y2, y1, y1))
	if (!is.null(xscale)) text((x1 + x2)/2, y1, xscale, pos=1)
	if (!is.null(yscale)) text(x1, (y1 + y2)/2, yscale, pos=4)
	# text(125, -69, "2 mV", pos=4)
}


scalebar2 <- function(dx, dy, xscale=NULL, yscale=NULL, pos='bottomleft'){
	lims <- par('usr')
	dxx <- lims[2]-lims[1]
	dyy <- lims[4]-lims[3]
	
	if (pos=='bottomleft'){
		x1 <- lims[1] + dxx/20
		x2 <- x1 + dx
		y1 <- lims[3] + dyy/20 + dy
		y2 <- y1 + dy
	}

	if (pos=='bottomright'){
		x1 <- lims[2] - dxx/20 - dx
		x2 <- x1 + dx
		y1 <- lims[3] + dyy/20 + dy
		y2 <- y1 + dy
	}

	if (pos=='topleft'){
		x1 <- lims[1] + dxx/20
		x2 <- x1 + dx
		y1 <- lims[4] - dyy/20 - dy
		y2 <- y1 + dy
	}

	if (pos=='topright'){
		x1 <- lims[2] - dxx/20 - dx
		x2 <- x1 + dx
		y1 <- lims[4] - dyy/20 - dy
		y2 <- y1 + dy
	}
	
	lines(c(x1, x1, x2), c(y2, y1, y1))
	if (!is.null(xscale)) text((x1 + x2)/2, y1, xscale, pos=1)
	if (!is.null(yscale)) text(x1, (y1 + y2)/2, yscale, pos=4)
	# text(125, -69, "2 mV", pos=4)
}


##########################################
## fitting data with a cosine function


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


####################################
### defining a colormap from grey - yellw - red

fr <- function(x){
	redvalue <- 0.9 + 1.41 * x - 1.81 * x^2
	redvalue[redvalue < 0] <- 0
	redvalue[redvalue > 1] <- 1
	redvalue
}

fg <- function(x){
	greenvalue <- 0.9 + 0.96 * x - 3.21 * x^2
	greenvalue[greenvalue < 0] <- 0
	greenvalue[greenvalue > 1] <- 1
	greenvalue
}

fb <- function(x){
	bluevalue <- 0.9 - 1.96 * x + 0.96 * x^2
	bluevalue[bluevalue < 0] <- 0
	bluevalue[bluevalue > 1] <- 1
	bluevalue
}

gyr <- function(N){
	x <- seq(0, N)/N
	x <- x^0.75
	cols <- rgb(fr(x), fg(x), fb(x), max=1)	
	cols
}

