### this is a script that calculates the mean decoding error as a function of the position shifted in space or in time
### its output is similar to Fig 2d, Fig 3f and Figure 4–Figure supplement 2.


infile <- paste('spikes_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, sampler, '.RData', sep='')	
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
## we similarly filter the position in the case of real data
txyvd <- posfilter(cbind(motiondata$T.axis, t(motiondata$position$x*100)), 0.5, F)
dt_pos <- diff(motiondata$T.axis[1:2])
txyvd <- rbind(txyvd, c(tail(txyvd, 1)[1] + dt_pos, tail(txyvd, 1)[2:3], 0, 0))
txy <- txyvd[,1:3]

start_br <- (min(txy[,2:3]) %/% 5) * 5
end_br <- (max(txy[,2:3]) %/% 5 +  1) * 5
xy_brs <- seq(start_br, end_br, by=5)
i_run <- seq(1, nrow(txy))

pmap <- estimate_posmap(txy, i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F)
istart <- min(which(attr(pmap, 'xcenter')>0))
iend <- max(which(attr(pmap, 'xcenter')<200))

ncells <- dim(data$rates)[1]
rmaps <- array(NA, dim=c(ncells, ncol(pmap), ncol(pmap)))

mapsname <- paste('3maps_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, '_', jitter, sampler, '.RData', sep='')	

if (mapsname %in% list.files()) {
	print('loading ratemaps...')
	load(file=mapsname)
	ratemaps <- maps$ratemaps
	# sharpmaps2 <- maps$sharpmaps[,2,,]
} else{
	par(mfcol=c(3,5))
	for (ii in 1:ncells){
		spi <- spt[spt[,2] == ii,1]
		ratemap <- estimate_ratemap(spi, txy, i_data=i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=pmap, graphics=0)
		rmaps[ii,,] <- ratemap
		
	}
	rmaps[rmaps<0.1] <- 0.1
	ratemaps <- rmaps[,istart:iend,istart:iend]
}

posmap <- pmap[istart:iend, istart:iend]
attr(posmap, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]
attr(posmap, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]

attr(ratemaps, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]
attr(ratemaps, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]

dir.create('./SimFigs/Fig3_FS3/', showW=F)
PFfilename <- paste('./SimFigs/Fig3_FS3/PFstats_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, evensp, '_', jitter, '.pdf', sep='')	
pdf(file=PFfilename, 10, 5, useD=F)
FLB <- plot_rate_stats(ratemaps, 1:N.cells, 1:N.cells, c(9, 10, 18, 19), posmap=posmap)
dev.off()

attr(ratemaps, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]/100
attr(ratemaps, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]/100		
attr(posmap, 'xcenters') <- attr(pmap, 'xcenter')[istart:iend]/100
attr(posmap, 'ycenters') <- attr(pmap, 'xcenter')[istart:iend]/100		



###################################################################
### MAP decoding the spikes in the three different window of the theta
###################################################################

outfile <- paste('decode_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, evensp, '_', jitter, sampler, '.RData', sep='')		

if (outfile %in% list.files()) {
	print('loading MAP decoded theta sequences...')
	load(file=outfile)
	postmean <- decode.theta$postmean
	postvar <- decode.theta$postvar
	nsp.theta <- decode.theta$nsp	
} else {
	print('MAP decoding theta sequences...')
	postmean <- array(NA, dim=c(2, 3, N.theta))
	postvar <- array(NA, dim=c(2, 3, N.theta)) 
	nsp.theta <- array(NA, dim=c(3, N.theta)) # spikes

	rasters <- array(0, dim=c(3, N.theta, N.cells)) # for posterior decoding
	bin_widths <- array(0, dim=c(3, N.theta))
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
			if (evenspikes&(nrow(sp.theta)>2)){
				L <- nrow(sp.theta)
				sp_phase <- sort(c(rep(1:3, floor(L / 3)), sample(1:3, 3))[1:L])
				L1 <- sum(sp_phase == 1)
				L2 <- sum(sp_phase == 2)
				L3 <- sum(sp_phase == 3)
				segbounds <- c(0, (sp.theta[L1,1] + sp.theta[L1+1,1])/2, (sp.theta[L1+L2,1] + sp.theta[L1+L2+1,1])/2, t2 - t1) ## uniform segments in spike counts
			} else {
				segbounds <- c(0, 1/3, 2/3, 1) * (t2 - t1) ## uniform segments in time		
			}
			bin_widths[, i.theta] <- diff(segbounds)
			for (i.seg in 1:3){ # past, present, future segments
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
	
	postmean[,1,] <- t(post1$mean)
	postmean[,2,] <- t(post2$mean)
	postmean[,3,] <- t(post3$mean)

	postvar[,1,] <- t(post1$var)
	postvar[,2,] <- t(post2$var)
	postvar[,3,] <- t(post3$var)
		
	decode.theta <- list(postmean=postmean, postvar=postvar, nsp=nsp.theta)
	save(decode.theta, file=outfile)
}


#########################################################
## rmaps based on early phase spikes with temporally shifted position data
#########################################################

mapsname <- paste('3maps_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, '_', jitter, sampler, '.RData', sep='')

if (mapsname %in% list.files()) {
	print('loading early, mid and late-based maps...')
	load(file=mapsname)
	ratemaps <- maps$ratemaps
	ratemaps_early <- maps$ratemaps_early
} else{
	## sorting spikes with respect to theta cycles
	print('calculating ratemaps from early spikes...')
	thetas <- list(cbind(theta.starts[1:N.theta], theta.starts[2:(N.theta+1)]))
	chunks <- get_theta_seqs(sp=spt[,1:2], txy=txyvd,  i_data=i_run, cellids=1:ncells, thetas=thetas)
	
	rast_early <- chunks$past[[1]]$rast
	txy_early <- chunks$past[[1]]$pos
	pmap_early <- estimate_posmap_xyw(txy_early[,c(2,3,5)], xbrs=xy_brs, ybrs=xy_brs, sigma=10, graphics=T, prior=1/100)
	rmaps_early <- array(NA, dim=c(ncells, nrow(pmap_early), ncol(pmap_early)))

	rast_late <- chunks$future[[1]]$rast
	txy_late <- chunks$future[[1]]$pos
	pmap_late <- estimate_posmap_xyw(txy_late[,c(2,3,5)], xbrs=xy_brs, ybrs=xy_brs, sigma=10, graphics=T, prior=1/100)
	rmaps_late <- array(NA, dim=c(ncells, nrow(pmap_late), ncol(pmap_late)))

	rast_mid <- chunks$present[[1]]$rast
	txy_mid <- chunks$present[[1]]$pos
	pmap_mid <- estimate_posmap_xyw(txy_mid[,c(2,3,5)], xbrs=xy_brs, ybrs=xy_brs, sigma=10, graphics=T, prior=1/100)
	rmaps_mid <- array(NA, dim=c(ncells, nrow(pmap_mid), ncol(pmap_mid)))
		
	### we need the decode.theta - decoded position based on the ratemaps to calculate this
	savefile <- F
	plot_enc_err <- F
	
	dir.create('./SimFigs/Fig3/', showW=F)
	curve_filename <- paste('./SimFigs/Fig3/decoding_shift_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, evensp, '_', jitter, '.pdf', sep='')	
	pdf(file= curve_filename, 8, 5, useD=F)
	source('./decoding_curves.R')
	dev.off()


	tau_shift <- (which.min(m.err[1,]) - 11)/10
	xx_fut <- approx(t.theta, x.theta*100, xout=t.theta + tau_shift)$y 
	yy_fut <- approx(t.theta, y.theta*100, xout=t.theta + tau_shift)$y 
	dd_fut <- approx(t.theta, d.theta, xout=t.theta + tau_shift)$y
	xyd <- cbind(xx_fut, yy_fut, dd_fut)
	
	for (ii in 1:ncells){
		counts <- rast_early[,ii]
		ratemap <- estimate_ratemap_xyw(counts, xyd, xbrs=xy_brs, ybrs=xy_brs, sigma=10, posmap=pmap_early, graphics=0)
		rmaps_early[ii,,] <- ratemap	
	}	
	ratemaps_early <- rmaps_early[,istart:iend,istart:iend]
	
	tau_shift <- (which.min(m.err[2,]) - 11)/10
	xx_fut <- approx(t.theta, x.theta*100, xout=t.theta + tau_shift)$y 
	yy_fut <- approx(t.theta, y.theta*100, xout=t.theta + tau_shift)$y 
	dd_fut <- approx(t.theta, d.theta, xout=t.theta + tau_shift)$y
	xyd <- cbind(xx_fut, yy_fut, dd_fut)
	
	for (ii in 1:ncells){
		counts <- rast_mid[,ii]
		ratemap <- estimate_ratemap_xyw(counts, xyd, xbrs=xy_brs, ybrs=xy_brs, sigma=10, posmap=pmap_mid, graphics=0)
		rmaps_mid[ii,,] <- ratemap	
	}	

	tau_shift <- (which.min(m.err[3,]) - 11)/10
	xx_fut <- approx(t.theta, x.theta*100, xout=t.theta + tau_shift)$y 
	yy_fut <- approx(t.theta, y.theta*100, xout=t.theta + tau_shift)$y 
	dd_fut <- approx(t.theta, d.theta, xout=t.theta + tau_shift)$y
	xyd <- cbind(xx_fut, yy_fut, dd_fut)
	
	for (ii in 1:ncells){
		counts <- rast_late[,ii]
		ratemap <- estimate_ratemap_xyw(counts, xyd, xbrs=xy_brs, ybrs=xy_brs, sigma=10, posmap=pmap_late, graphics=0)
		rmaps_late[ii,,] <- ratemap	
	}
	
	maps <- list(ratemaps=ratemaps, ratemaps_early=ratemaps_early)			

	maps$ratemaps_late <- rmaps_late[,istart:iend,istart:iend]
	maps$ratemaps_mid <- rmaps_mid[,istart:iend,istart:iend]
	save(maps, file=mapsname)
}


##############################################
## decoding shifts - only with MAP decoding
##############################################
## we need to load the correct decode.theta

decfile <- paste('decode_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, evensp, '_', jitter, sampler, '.RData', sep='')		
load(file=decfile)

source('decoding_curves.R')

