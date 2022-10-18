#  DDC decoding using standard ratemaps
dir.create('./SimFigs/Fig5/', showW=F)

## decoding the spikes 
outfile <- paste('decode_estR_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, '_DDC.RData', sep='')		
if (outfile %in% list.files('./SimFigs/Fig5/')) {
	ddcfile <- paste('./SimFigs/Fig5/', outfile, sep='')
	load(file=ddcfile)
	opars <- decode_DDC$opars
} else {

	thetas <- list(cbind(theta.starts[1:N.theta], theta.starts[2:(N.theta+1)]))
	chunks <- get_theta_seqs(sp=spt[,1:2], txy=txyvd,  i_data=i_run, cellids=1:ncells, thetas=thetas)
	
	rast_early <- chunks$past[[1]]$rast
	rast_mid <- chunks$present[[1]]$rast
	rast_late <- chunks$future[[1]]$rast

	txy_early <- chunks$past[[1]]$pos
	txy_mid <- chunks$present[[1]]$pos
	txy_late <- chunks$future[[1]]$pos

	bincenters <- attr(posmap, 'xcenters')
	Nbins <- length(bincenters)
	
	g.x <- matrix(rep(bincenters, Nbins), Nbins) * 100
	g.y <- matrix(rep(bincenters, each=Nbins), Nbins) * 100
	
	g.x.vec <- as.vector(g.x)
	g.y.vec <- as.vector(g.y)
	
	rates.vec <- matrix(0, Nbins*Nbins, ncells)

	for (i in 1:ncells){
		rates.vec[,i] <- as.vector(ratemaps[i,,])	
	}
	rmin <- 0.01 # min firing rate
	rates.vec[rates.vec < rmin] <- rmin
	
	# i.cycle <- 1
	# t_theta <- txy_early[i.cycle,5]
	# sp <- rast_early[i.cycle,]
	
	initpars <- c(mu1=100, mu2=100, sigma=5)
	# DDC_likelihood(initpars, sp, rates.vec, g.x.vec, g.y.vec, t_theta)
	
	opars <- array(NA, dim=c(3, 4, N.theta), dimnames=list(c('early', 'mid', 'late'), c('mean1', 'mean2', 'sd', 'nsp'), NULL))
	
	for (i.cycle in 1:N.theta){
		sp <- rast_early[i.cycle,]
		t_theta <- txy_early[i.cycle,5]
		opars[1,, i.cycle] <- c(optim(initpars, DDC_likelihood, gr=NULL, sp=sp, rates_vec=rates.vec, g.x.vec=g.x.vec, g.y.vec=g.y.vec, deltaT=t_theta, method='L-BFGS-B', lower=c(0, 0, 3), upper=c(200, 200, 200), control=list(fnscale=-1))$par, sum(sp))

		sp <- rast_mid[i.cycle,]
		t_theta <- txy_mid[i.cycle,5]
		opars[2,, i.cycle] <- c(optim(initpars, DDC_likelihood, gr=NULL, sp=sp, rates_vec=rates.vec, g.x.vec=g.x.vec, g.y.vec=g.y.vec, deltaT=t_theta, method='L-BFGS-B', lower=c(0, 0, 3), upper=c(200, 200, 200), control=list(fnscale=-1))$par, sum(sp))

		sp <- rast_late[i.cycle,]
		t_theta <- txy_late[i.cycle,5]
		opars[3,,i.cycle] <- c(optim(initpars, DDC_likelihood, gr=NULL, sp=sp, rates_vec=rates.vec, g.x.vec=g.x.vec, g.y.vec=g.y.vec, deltaT=t_theta, method='L-BFGS-B', lower=c(0, 0, 3), upper=c(200, 200, 200), control=list(fnscale=-1))$par, sum(sp))
		if (i.cycle%%100 == 0) cat(i.cycle, '\n')
	}
	decode_DDC <- list(opars=opars)
	ddcfile <- paste('./SimFigs/Fig5/', outfile, sep='')
	save(decode_DDC, file=ddcfile)
}
