########################################################################
## combine LFP, position and spike data into a single, interactive plot

#######################################################################		
## preparing the ratemaps - we need them in a vectorial form...
print('preparing data for DDC decoding')


mapsname <- paste('3maps_rat', i_rat, '_day', i_day, '.RData', sep='')
if (mapsname %in% list.files()) {
	print('loading ratemaps')
	load(file=mapsname)
	ratemaps <- maps$ratemaps
	ratemaps_eml <- maps$ratemaps_eml
} else {
	print('mapfile not found! Run MAP_decode.R first!')
}

outfile <- paste('R', i_rat, 'D', i_day, '_decoded_theta_seqs.RData', sep='')
if (outfile %in% list.files()) {
	print('loading decoded position')
	load(file=outfile)
	chunks <- dec_theta$chunks
	estimpos_early <- dec_theta$est_early
	estimpos_mid <- dec_theta$est_mid
	estimpos_late <- dec_theta$est_late
} else {
	print('decoded locations not found! Run MAP_decode.R first!')
}


bincenters <- attr(posmap, 'xcenters')
Nbins <- length(bincenters)
g.x <- matrix(rep(bincenters, Nbins), Nbins)
g.y <- matrix(rep(bincenters, each=Nbins), Nbins)

g.x.vec <- as.vector(g.x)
g.y.vec <- as.vector(g.y)

rates.vec <- matrix(0, Nbins*Nbins, ncells)
for (i in 1:ncells){
	rates.vec[,i] <- as.vector(ratemaps[i,,])	
}

rmin <- 0.1 # min firing rate
rates.vec[rates.vec < rmin] <- rmin

initpars <- c(mu1=100, mu2=100, sigma=100)

#######################################################################		
## preparing the spike rasters ...

rast_early <- chunks$past[[1]]$rast
txy_early <- chunks$past[[1]]$pos

rast_mid <- chunks$present[[1]]$rast
txy_mid <- chunks$present[[1]]$pos

rast_late <- chunks$future[[1]]$rast
txy_late <- chunks$future[[1]]$pos

for (i_chunk in 2:length(chunks$past)){
	rast_early <- rbind(rast_early, chunks$past[[i_chunk]]$rast)
	txy_early <- rbind(txy_early, chunks$past[[i_chunk]]$pos)

	rast_mid <- rbind(rast_mid, chunks$present[[i_chunk]]$rast)
	txy_mid <- rbind(txy_mid, chunks$present[[i_chunk]]$pos)

	rast_late <- rbind(rast_late, chunks$future[[i_chunk]]$rast)
	txy_late <- rbind(txy_late, chunks$future[[i_chunk]]$pos)
}
	
#######################################################################		
## decoding the spike rasters ...
print('starting DDC decoding')

NN <- nrow(rast_early)

initpars <- c(mu1=100, mu2=100, sigma=5)
i.cycle <- 1
t_theta <- txy_early[i.cycle,5]
sp <- rast_early[i.cycle,]
DDC_likelihood(initpars, sp, rates.vec, g.x.vec, g.y.vec, t_theta)

oparfile <- paste('R', i_rat, 'D', i_day, '_DDC_MLpars_th3.RData', sep='')
if (oparfile %in% list.files()){
	load(oparfile)
	opars <- MLpars$opars
	truepos <- MLpars$truepos
} else {
	opars <- array(NA, dim=c(3, 4, NN), dimnames=list(c('early', 'mid', 'late'), c('mean1', 'mean2', 'sd', 'nsp'), NULL))
	truepos <- array(NA, dim=c(3, 4, NN), dimnames=list(c('early', 'mid', 'late'), c('x', 'y', 'dir', 'binwidth'), NULL))
	
	for (i.cycle in 15818:NN){
		sp <- rast_early[i.cycle,]
		opars[1,4,i.cycle] <- sum(sp, na.rm=T)
	
		t_theta <- txy_early[i.cycle,5]
		opars[1,1:3,i.cycle] <- optim(initpars, DDC_likelihood, gr=NULL, sp=sp, rates_vec=rates.vec, g.x.vec=g.x.vec, g.y.vec=g.y.vec, deltaT=t_theta, method='L-BFGS-B', lower=c(0, 0, 3), upper=c(200, 200, 200), control=list(fnscale=-1))$par
		truepos[1,,i.cycle] <- txy_early[i.cycle,2:5]
	
		sp <- rast_mid[i.cycle,]
		opars[2,4,i.cycle] <- sum(sp, na.rm=T)
	
		t_theta <- txy_mid[i.cycle,5]
		opars[2,1:3,i.cycle] <- optim(initpars, DDC_likelihood, gr=NULL, sp=sp, rates_vec=rates.vec, g.x.vec=g.x.vec, g.y.vec=g.y.vec, deltaT=t_theta, method='L-BFGS-B', lower=c(0, 0, 3), upper=c(200, 200, 200), control=list(fnscale=-1))$par
		truepos[2,,i.cycle] <- txy_mid[i.cycle,2:5]
	
		sp <- rast_late[i.cycle,]
		opars[3,4,i.cycle] <- sum(sp, na.rm=T)
	
		t_theta <- txy_late[i.cycle,5]
		opars[3,1:3,i.cycle] <- optim(initpars, DDC_likelihood, gr=NULL, sp=sp, rates_vec=rates.vec, g.x.vec=g.x.vec, g.y.vec=g.y.vec, deltaT=t_theta, method='L-BFGS-B', lower=c(0, 0, 3), upper=c(200, 200, 200), control=list(fnscale=-1))$par
		truepos[3,,i.cycle] <- txy_late[i.cycle,2:5]
		if (i.cycle %% 100 == 0) cat(i.cycle, ' ')
	}
	
	MLpars <- list(opars=opars, truepos=truepos)
	save(MLpars, file=oparfile)
}
###########################################################################
## plot the data

th.nsp <- 15
mean(opars[,4,])
dir.create('./DataFigures/Fig5/', showW=F)

for (th.nsp in c(0, median(opars[,4,]))){
	th_SD <- 50
	
	ii.gr1 <- opars[1,4,] > th.nsp
	ii.gr2 <- opars[2,4,] > th.nsp
	ii.gr3 <- opars[3,4,] > th.nsp
	
	CDF_SD <- cbind(cumsum(hist(opars[1,3,ii.gr1], breaks=3:201-0.5, plot=F)$density), cumsum(hist(opars[2,3,ii.gr2], breaks=3:201-0.5, plot=F)$density), cumsum(hist(opars[3,3,ii.gr3], breaks=3:201-0.5, plot=F)$density))
	
	ms_o <- rep(NA, 3)
	ss_o <- rep(NA, 3)

	ms_o_th <- rep(NA, 3)
	ss_o_th <- rep(NA, 3)	
	
	i1 <- opars[1,3,] < th_SD
	i2 <- opars[2,3,] < th_SD
	i3 <- opars[3,3,] < th_SD
	ks.test(opars[1,3,i1&ii.gr1], opars[2,3, i2&ii.gr2], alternative='greater')
	Pa1 <- round(log(ks.test(opars[1,3,i1& ii.gr1], opars[3,3,i3&ii.gr3], alternative='greater')$p, 10))
	Pa2 <- round(log(ks.test(opars[2,3,i2& ii.gr2], opars[3,3,i3& ii.gr3], alternative='greater')$p, 10))
	Pa <- min(c(Pa1, Pa2))
		
	{
		ms_o[1] <- mean(opars[1,3, ii.gr1])
		ss_o[1] <- sd(opars[1,3, ii.gr1]) / sqrt(sum(ii.gr1))
		
		ms_o[2] <- mean(opars[2,3, ii.gr2])
		ss_o[2] <- sd(opars[2,3,ii.gr2]) / sqrt(sum(ii.gr2))
		
		ms_o[3] <- mean(opars[3,3,ii.gr3])
		ss_o[3] <- sd(opars[3,3,ii.gr3]) / sqrt(sum(ii.gr3))
	}
	
	{
		ms_o_th[1] <- mean(opars[1,3,i1& ii.gr1])
		ss_o_th[1] <- sd(opars[1,3,i1& ii.gr1]) / sqrt(sum(i1& ii.gr1))
		
		ms_o_th[2] <- mean(opars[2,3,i2& ii.gr2])
		ss_o_th[2] <- sd(opars[2,3,i2& ii.gr2]) / sqrt(sum(i2& ii.gr2))
		
		ms_o_th[3] <- mean(opars[3,3,i3& ii.gr3])
		ss_o_th[3] <- sd(opars[3,3,i3& ii.gr3]) / sqrt(sum(i3& ii.gr3))
	}
	
	DDCfile <- paste('DataFigures/Fig5/R', i_rat, 'D', i_day, 'DDC_decSD_smap1_nsp', round(th.nsp,1), '_th3.pdf', sep='')		
	pdf(file=DDCfile, 5, 2.5, useD=F)
	
	par(mfcol=c(1,2))
	par(mar=c(4,4,1,1))
	
	matplot(3:200, CDF_SD, t='l', ylim=c(0.8, 1), lty=1, col=viridis(5, option='B')[2:4], lwd=c(4, 2,2 ), xlim=c(0, 50), xlab='SD (cm)', ylab='cumulative probability', main=paste('all', Pa))
	legend('bottomright', leg=c('early', 'mid', 'late'), lty=1, lwd=c(4,2,2), col=viridis(5, option='B')[2:4], bty='n')
	
	
	plot(0:2*4, ms_o, xlim=c(-1, 12), axes=F, ylim=c(3, 12), cex=2, t='o', ylab='standard deviation', xlab='', pch=21, col=viridis(5, option='B')[2:4], bg=grey(1))
	lines(0:2*4, ms_o_th, t='o', pch=21, bg=viridis(5, option='B')[2:4], cex=2)
	
	plotCI(0:2*4, ms_o, ss_o, gap=0, add=T, t='l')
	plotCI(0:2*4, ms_o_th, ss_o_th, gap=0, add=T, t='l')
	
	axis(1, 0:2*4+0.5, c('early', 'mid', 'late'), tick=F)
	axis(2, c(3, 6, 9, 12), las=2)
	legend('topleft', leg=c('ratemap, all',  'SD<50'), pch=21, pt.cex=c(1.5), col=c(viridis(5, option='B')[3], 1), pt.bg=c(grey(1), viridis(5, option='B')[3]), bty='n')
	dev.off()
}