## loading the spikes 
infile <- paste('spikes_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, '.RData', sep='')
load(file=infile)
txy_mu <- data$xy.trajectories # the posterior mean (100 Hz)
txy_var <- data$var.trajectories # the posterior variance (100 Hz)

## loading the motion data
infile <- paste('motiondata_Tmax', Tmax, reg, '.RData', sep='')
load(file=infile)


### spikes	
spt <- data$spt
mean.rate <- nrow(spt) / N.cells / diff(range(spt[,1]))
cat('average firing rate: ', mean.rate, 'Hz, number of spikes within theta: ', mean.rate * N.cells * 0.1) # 2.4 Hz / cell - 24 spikes/theta
theta.starts.orig <- data$theta.starts
N.theta <- length(theta.starts.orig)-1


jitterfile <- paste('jitter_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, evensp, '_', jitter, '.RData', sep='')
if (jitterfile %in% list.files()) {
	load(jitterfile)
} else {
	jitter_ms <- jitter / 1000/2
	randmin <- -jitter_ms
	randmax <- jitter_ms	
	theta.starts <- theta.starts.orig + runif(N.theta + 1, randmin, randmax)
	save(theta.starts, file=jitterfile)
}

## filtering the position data and adding velocity and direction to position info
## we do filter the position in the case of real data
txyvd <- posfilter(cbind(motiondata$T.axis, t(motiondata$position$x*100)), 0.5, F)
dt_pos <- diff(motiondata$T.axis[1:2])
txyvd <- rbind(txyvd, c(tail(txyvd, 1)[1] + dt_pos, tail(txyvd, 1)[2:3], 0, 0))
txy <- txyvd[,1:3]

start_br <- (min(txy[,2:3]) %/% 5) * 5
end_br <- (max(txy[,2:3]) %/% 5 +  1) * 5
xy_brs <- seq(start_br, end_br, by=5)
i_run <- seq(1, nrow(txy))

pmap <- estimate_posmap(txy, i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F)

ncells <- dim(data$rates)[1]
rmaps <- array(NA, dim=c(ncells, ncol(pmap), ncol(pmap)))
nspikes <- rep(NA, ncells)

par(mfcol=c(3,5))
for (ii in 1:ncells){
	spi <- spt[spt[,2] == ii,1]
	ratemap <- estimate_ratemap(spi, txy, i_data=i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=pmap, graphics=0)
	rmaps[ii,,] <- ratemap
	nspikes[ii] <- attr(ratemap, 'nsp')
	
}
rmaps[rmaps<0.1] <- 0.1

istart <- min(which(attr(pmap, 'xcenter')>0))
iend <- max(which(attr(pmap, 'xcenter')<200))
ratemaps <- rmaps[,istart:iend,istart:iend]
attr(ratemaps, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]
attr(ratemaps, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]

posmap <- pmap[istart:iend, istart:iend]
attr(posmap, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]
attr(posmap, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]


placestatfile <- paste('SimFigs/00PFstats/PFstats_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, evensp, '_', jitter, '.png', sep='')
if (savefile) png(file=placestatfile, 1200, 600, pointsize=18)
FLB <- plot_rate_stats(ratemaps, 1:N.cells, 1:N.cells, c(9, 10, 18, 19), posmap=posmap)
if (savefile) dev.off()

attr(ratemaps, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]/100
attr(ratemaps, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]/100		
attr(posmap, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]/100
attr(posmap, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]/100		

###################################################################
### MAP decoding the spikes in the 6 different window of the theta

## decoding the spikes 
outfile <- paste('decode6Eq_', code, reg, evensp, '_', jitter, '.RData', sep='')		

if (outfile %in% list.files()) {
	load(file=outfile)
	postmean <- decode.theta$postmean
	postvar <- decode.theta$postvar
	nsp.theta <- decode.theta$nsp	
} else {
	postmean <- array(NA, dim=c(2, 6, N.theta))
	postvar <- array(NA, dim=c(2, 6, N.theta)) 
	nsp.theta <- array(NA, dim=c(6, N.theta)) # spikes

	rasters <- array(0, dim=c(6, N.theta, N.cells)) # for posterior decoding
	bin_widths <- array(0, dim=c(6, N.theta))
	n2 <- 0
	t1 <- 0
	for (i.theta in 1:N.theta){
		t2 <- theta.starts[i.theta+1]
		t1 <- theta.starts[i.theta]
		n1 <- n2 + 1
		nn <- min(n1+5000, nrow(spt))
		nnn <- max(which(spt[n1:nn,1] < t2))
		if (nnn >= 5000) warning('too many spikes!')
		n2 <- n1 + nnn - 1
		if (n2 > n1) {
			sp.theta <- spt[n1:n2,1:2]
			sp.theta[,1] <- sp.theta[,1] - t1
			segbounds <- 0:6/6 * (t2 - t1) ## uniform segments in time		
			bin_widths[, i.theta] <- diff(segbounds)
			for (i.seg in 1:6){ # past, present, future segments
				ii.seg <- which((sp.theta[,1] >= segbounds[i.seg]) & (sp.theta[,1] < segbounds[i.seg+1]))
				sp <- sp.theta[ii.seg,2]
				ss <- spt2spc(sp, N.cells)
				rasters[i.seg, i.theta,] <- ss # for posterior decoding
				nsp.theta[i.seg, i.theta] <- length(sp)
			}
		}
		if (i.theta %% 1000 == 0) cat(i.theta, ' ')
	}

	post1 <- decode_ratemaps(rasters[1,,], bin_widths[1,], ratemaps, posmap)
	post2 <- decode_ratemaps(rasters[2,,], bin_widths[2,], ratemaps, posmap)
	post3 <- decode_ratemaps(rasters[3,,], bin_widths[3,], ratemaps, posmap)
	post4 <- decode_ratemaps(rasters[4,,], bin_widths[4,], ratemaps, posmap)
	post5 <- decode_ratemaps(rasters[5,,], bin_widths[5,], ratemaps, posmap)
	post6 <- decode_ratemaps(rasters[6,,], bin_widths[6,], ratemaps, posmap)
	
	postmean[,1,] <- t(post1$mean)
	postmean[,2,] <- t(post2$mean)
	postmean[,3,] <- t(post3$mean)
	postmean[,4,] <- t(post4$mean)
	postmean[,5,] <- t(post5$mean)
	postmean[,6,] <- t(post6$mean)

	postvar[,1,] <- t(post1$var)
	postvar[,2,] <- t(post2$var)
	postvar[,3,] <- t(post3$var)
	postvar[,4,] <- t(post4$var)
	postvar[,5,] <- t(post5$var)
	postvar[,6,] <- t(post6$var)
		
	decode.theta <- list(postmean=postmean, postvar=postvar, nsp=nsp.theta)
	save(decode.theta, file=outfile)
}


##############################################################
## compare the decoded and the reference position

t.theta <- data$theta.starts[1:N.theta]
x.theta <- approx(motiondata$T.axis, txyvd[,2], t.theta)$y / 100
y.theta <- approx(motiondata$T.axis, txyvd[,3], t.theta)$y / 100
v.theta <- approx(motiondata$T.axis, txyvd[,4], t.theta)$y
d.theta <- approx(motiondata$T.axis, txyvd[,5], t.theta)$y
xy.theta <- cbind(x.theta, y.theta)

relative_xyerr <- postmean
for (i.t in 1:6){
	dxi <- postmean[1,i.t,] - x.theta
	dyi <- postmean[2,i.t,] - y.theta
	rel_direction <- atan2(dyi, dxi) - d.theta
	err <- sqrt(dxi^2 + dyi^2)
	relative_xyerr[1, i.t, ] <- err * sin(rel_direction)
	relative_xyerr[2, i.t, ] <- err * cos(rel_direction)
}

xvar <- rowMeans(postvar[2,,]) * 10000
yvar <- rowMeans(postvar[1,,]) * 10000


library(colormap)
cols <- colormap(colormaps$phase, 7)[1:6]
outfile <- paste('SimFigs/07phase6/decode6Eq_', code, reg, evensp, '_', jitter, '.pdf', sep='')		
if (savefile) pdf(file=outfile, 6, 6, useD=F)
par(mfcol=c(2,1))
par(mar=c(1,4,1,1))
plot(1:12, rep(rowMeans(relative_xyerr[2,,])*100, 2), pch=21, bg=cols, cex=2, t='o', xlab='theta phase (deg)', ylab='mean error (cm)', axes=F, xlim=c(0, 13))
points(1:12, rep(rowMeans(relative_xyerr[1,,])*100, 2), pch=22, bg=cols, cex=2, t='o')
axis(2, las=2)


par(mar=c(4,4,1,1))
plot(1:12, rep(xvar, 2), pch=21, bg=cols, cex=1.5, t='o', xlab='theta phase (deg)', ylab='posterior variance', axes=F, xlim=c(0, 13), ylim=c(600, 1100))# ylim=range(c(xvar, yvar)))
points(1:12, rep(yvar, 2), pch=22, bg=cols, cex=1.5, t='o')
axis(2, las=2)
axis(1, c(0, 3, 6, 9, 12) + 0.5, c(0, 180, 360, 180, 360))
if (savefile) dev.off()

