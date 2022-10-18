

require(abind)

### combine position, decoded position and number of spikes across theta cycles
txyd_theta <- abind(chunks$past[[1]]$pos, chunks$present[[1]]$pos, chunks$future[[1]]$pos, along=-1)
mu_theta <- abind(estimpos_early[[1]]$estim_pos, estimpos_mid[[1]]$estim_pos, estimpos_late[[1]]$estim_pos, along=-1)
nsp <- rowSums(chunks$present[[1]]$rast, na.rm=T)


for (i_chunk in 2:length(chunks$past)){
	if (length(which(!is.na(chunks$past[[i_chunk]]$pos[,1]))) > 0){
		txyd_theta_i <- abind(chunks$past[[i_chunk]]$pos, chunks$present[[i_chunk]]$pos, chunks$future[[i_chunk]]$pos, along=-1)
		mu_theta_i <- abind(estimpos_early[[i_chunk]]$estim_pos, estimpos_mid[[i_chunk]]$estim_pos, estimpos_late[[i_chunk]]$estim_pos, along=-1)
		txyd_theta <- abind(txyd_theta, txyd_theta_i, along=2)
		mu_theta <- abind(mu_theta, mu_theta_i, along=2)
		nsp <- c(nsp, rowSums(chunks$present[[i_chunk]]$rast, na.rm=T))
	}
}

N.theta <- dim(txyd_theta)[2]

quantile_curvature <- 0
th_curvature <- 0

## calculate the errors wrt the temporally or spatially shifted position
errors <- array(NA, dim=c(3, 31, N.theta))
errors_fw <- array(NA, dim=c(3, 31, N.theta)) # forward prediction

i.eml <- 1; i.tau <- 1

for (i.eml in 1:3){
	for (i.tau in 1:31){	
		## forward
		xx_fwd <- txyd_theta[i.eml,,2] + cos(txyd_theta[2,,4]) * 2 * (i.tau-11)
		yy_fwd <- txyd_theta[i.eml,,3] + sin(txyd_theta[2,,4]) * 2 * (i.tau-11)
		errors_fw[i.eml, i.tau, ] <- sqrt(rowSums((cbind(xx_fwd, yy_fwd) - mu_theta[i.eml,,])^2))
		
		## predict
		xx_fut <- approx(txyvd[,1], txyvd[,2], xout=txyd_theta[2,,1] + (i.tau-11)/10)$y
		yy_fut <- approx(txyvd[,1], txyvd[,3], xout=txyd_theta[2,,1] + (i.tau-11)/10)$y
		errors[i.eml, i.tau, ] <- sqrt(rowSums((cbind(xx_fut, yy_fut) - mu_theta[i.eml,,])^2))
	}
}

ths <- c(0, median(nsp), quantile(nsp, 7/8), quantile(nsp, 19/20))
m_enc_er <- array(NA, c(3, 4, 3), dimnames=list(c('mean', 'mean-outl', 'median'), paste('threshold', ths), c('start', 'mid', 'end')))
se_enc_er <- array(NA, c(3, 4, 3), dimnames=list(c('mean', 'mean-outl', 'median'), paste('threshold', ths), c('start', 'mid', 'end')))


tt <- -10:20

pdf(file=shiftfile, 8, 5, useD=F)
par(mfcol=c(1,1))
for (k in 4:1){
	nsp_th <- ths[k]
	
	m.err_fw <- matrix(NA, 3, 31)
	m.err <- matrix(NA, 3, 31)
	
	ii_nsp <- (nsp > nsp_th)
	for (i.eml in 1:3){
		m.err_fw[i.eml,] <- apply(errors_fw[i.eml,,ii_nsp], 1, median, na.rm=T) / 100# from cm to m^2
		m.err[i.eml,] <- apply(errors[i.eml,,ii_nsp], 1, median, na.rm=T) / 100
	}
	
	if (k ==4) cols <- viridis(10,option=3)[c(1,5,7)] else cols <- viridis(10,option=3, alpha=k/6)[c(1,5,7)]
	if (k==4){
		matplot(-10:20, t(m.err_fw), col=cols, pch=15:17, t='o', lty=1, ylim=c(0.05, 0.35), xlim=c(-10, 55), axes=F, xlab='', ylab='median error (m)')
		matplot(-10:20+35, t(m.err), col=cols, pch=15:17, t='o', lty=1, add=T)
		abline(h=apply(m.err, 1, min), col=cols, lty=3)
		legend('topleft', c('start', 'mid', 'end'), lty=1, pch=15:17, col=cols, bty='n')
		legend('topright', c(paste('th =', ths)), lty=1, col=viridis(1, option=3, alpha=c(1/6, 2/6, 3/6, 1)), bty='n')
		abline(v=c(0, 35))
		
		axis(1, seq(-10, 20, by=5)+35, seq(-1, 2, by=0.5))
		mtext('prediction (s)', 1, line=2,adj=0.8)
		
		axis(1, seq(-10, 20, by=5), seq(-20, 40, by=10))
		mtext('forward (cm)', 1, line=2,adj=0.2)
		# axis(2, 1:3/10, 1:3*10, las=2)
		axis(2, las=2)
	} else {
		matplot(-10:20, t(m.err_fw), col=cols, t='l', lty=1, add=T)
		matplot(-10:20+35, t(m.err), col=cols, t='l', lty=1, add=T)
	}
}

dev.off()

delay_eml <- tt[apply(m.err, 1, which.min)] / 10 # the delay of early, mid and late phase encoded position, in s
