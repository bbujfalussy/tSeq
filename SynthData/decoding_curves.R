if (jitter == 0){
	theta.starts <- data$theta.starts
	N.theta <- length(theta.starts)-1
}

# x and y position at theta start times
t.theta <- data$theta.starts[1:N.theta]
x.theta <- approx(txyvd[,1], txyvd[,2], t.theta)$y / 100
y.theta <- approx(txyvd[,1], txyvd[,3], t.theta)$y / 100
v.theta <- approx(txyvd[,1], txyvd[,4], t.theta)$y
d.theta <- approx(txyvd[,1], txyvd[,5], t.theta)$y
xy.theta <- cbind(x.theta, y.theta)


th_curvature <- 0
jj_curvature <- rep(T, N.theta)

ths <- c(0, median(nsp.theta), quantile(nsp.theta, 7/8), quantile(nsp.theta, 19/20))
errors <- array(NA, dim=c(3, 21, N.theta)) # temporal prediction
errors_fw <- array(NA, dim=c(3, 21, N.theta)) # forward prediction

for (i.ppf in 1:3){
	for (i.tau in 1:21){	
		## forward
		xx_fwd <- x.theta + cos(d.theta) * 0.02 * (i.tau-11)
		yy_fwd <- y.theta + sin(d.theta) * 0.02 * (i.tau-11)
		errors_fw[i.ppf, i.tau, ] <- sqrt(colSums((rbind(xx_fwd, yy_fwd) - decode.theta$postmean[,i.ppf,])^2))
		
		## predict
		xx_fut <- approx(t.theta, x.theta, xout=t.theta + (i.tau-11)/10)$y
		yy_fut <- approx(t.theta, y.theta, xout=t.theta + (i.tau-11)/10)$y
		errors[i.ppf, i.tau, ] <- sqrt(colSums((rbind(xx_fut, yy_fut) - decode.theta$postmean[,i.ppf,])^2))
	}
}


m_enc_er <- array(NA, c(3, 4, 3), dimnames=list(c('mean', 'mean-outl', 'median'), paste('threshold', ths), c('start', 'mid', 'end')))
se_enc_er <- array(NA, c(3, 4, 3), dimnames=list(c('mean', 'mean-outl', 'median'), paste('threshold', ths), c('start', 'mid', 'end')))

par(mfcol=c(1,1))
par(mar=c(4,4, 1, 1))

tt <- -10:10

TEE <- matrix(NA, 3, dim(errors)[3])
TEE_all <- list()
TEE_median <- list()

for (k in 4:1){
	nsp_th <- ths[k]
	
	m.err_fw <- matrix(NA, 3, 21)
	m.err <- matrix(NA, 3, 21)

	if (k ==4) cols <- viridis(10,option=3)[c(1,5,7)] else cols <- viridis(10,option=3, alpha=k/6)[c(1,5,7)]
	
	for (i.ppf in 1:3){
		if (i.ppf == 1) ii_nsp <- (decode.theta$nsp[i.ppf,] > nsp_th) & jj_curvature
		m.err_fw[i.ppf,] <- apply(errors_fw[i.ppf,,ii_nsp], 1, mean, na.rm=T)
		m.err[i.ppf,] <- apply(errors[i.ppf,,ii_nsp], 1, mean, na.rm=T)
		
		imin <- which.min(m.err[i.ppf,])
		q95 <- quantile(errors[i.ppf,imin,], 0.95, na.rm=T)
		ii_nsp95 <- errors[i.ppf,imin,] < q95
		if (k==1) {
			TEE[i.ppf,] <- (errors[i.ppf,imin,])^2
			TEE_all[[i.ppf]] <- (errors[i.ppf,imin,])^2
		}
		if (k==2) {
			TEE_median[[i.ppf]] <- (errors[i.ppf,imin,ii_nsp])^2
		}
		
		m_enc_er[1, k, i.ppf] <- mean((errors[i.ppf,imin,ii_nsp])^2, na.rm=T) 
		se_enc_er[1, k, i.ppf] <- sd((errors[i.ppf,imin,ii_nsp])^2, na.rm=T) / sqrt(length(ii_nsp))

		m_enc_er[2, k, i.ppf] <- mean((errors[i.ppf,imin,ii_nsp&ii_nsp95])^2, na.rm=T) # from cm to m^2
		se_enc_er[2, k, i.ppf] <- sd((errors[i.ppf,imin,ii_nsp&ii_nsp95])^2, na.rm=T) / sqrt(length(ii_nsp))

		m_enc_er[3, k, i.ppf] <- median((errors[i.ppf,imin,ii_nsp])^2, na.rm=T) # from cm to m^2
		se_enc_er[3, k, i.ppf] <- sd((errors[i.ppf,imin,ii_nsp])^2, na.rm=T) / sqrt(length(ii_nsp))
	}
	
	if (k==4){
		matplot(-10:10, t(m.err_fw), col=cols, pch=15:17, t='o', lty=1, ylim=c(0.05, 0.35), xlim=c(-10, 35), axes=F, xlab='', ylab='mean error (m)')
		matplot(-10:10+25, t(m.err), col=cols, pch=15:17, t='o', lty=1, add=T)
		abline(h=apply(m.err, 1, min), col=cols, lty=3)
		legend('topleft', c('past', 'present', 'future'), lty=1, pch=15:17, col=cols, bty='n')
		legend('topright', c(paste('spike threshold =', ths)), lty=1, col=viridis(1, option=3, alpha=c(1/6, 2/6, 3/6, 1)), bty='n')
		abline(v=c(0, 25))
		
		axis(1, seq(-10, 10, by=5)+25, seq(-1, 1, by=0.5))
		mtext('prediction (s)', 1, line=2,adj=0.8)
		
		axis(1, seq(-10, 10, by=5), seq(-20, 20, by=10))
		mtext('forward (m)', 1, line=2,adj=0.2)
		# axis(2, 1:3/10, 1:3*10, las=2)
		axis(2, las=2)
	} else {
		matplot(-10:10, t(m.err_fw), col=cols, t='l', lty=1, add=T)
		matplot(-10:10+25, t(m.err), col=cols, t='l', lty=1, add=T)
	}
}

delay_eml <- tt[apply(m.err, 1, which.min)] / 10 # the delay of early, mid and late phase encoded position, in s
