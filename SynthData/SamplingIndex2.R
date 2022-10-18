
## calculates cycle-to-cycle variability and the sampling index
#####################################################
## loading decoded position
#####################################################
se <- function(x, ...) sqrt(var(x, ...)/sum(!is.na(x)))
dir.create('./SimFigs/Fig6/', showW=F)


decfile <- paste('decode_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, evensp, '_', jitter, sampler, '.RData', sep='')		
load(file=decfile)

## we calculate the encoding error here...
quantile_curvature <- 0

sf <- savefile
savefile <- F
source('decoding_curves.R')
savefile <- sf


#####################################################
## calculating reference position
#####################################################

# xz coordinates from the encoding model - needed to check the calculations
t.theta <- data$theta.starts[1:N.theta]
x.theta <- approx(txyvd[,1], txyvd[,2], t.theta)$y / 100
y.theta <- approx(txyvd[,1], txyvd[,3], t.theta)$y / 100
xy_theta_true <- array(NA, dim=c(2, 3, N.theta)) # true position
for (i.tau in 1:3){
	delay <- delay_eml[i.tau]
	xy_theta_true[1,i.tau,] <- approx(t.theta, x.theta, xout=t.theta + delay)$y
	xy_theta_true[2,i.tau,] <- approx(t.theta, y.theta, xout=t.theta + delay)$y
}

xy_theta_repr <- array(NA, dim=c(2, 3, N.theta)) # represented position
mu_theta_post <- array(NA, dim=c(2, 3, N.theta)) # posterior mean - different from xy_theta_repr only in sampling!
var_theta_post <- array(NA, dim=c(2, 3, N.theta)) # posterior variance	

N.pos <- ncol(motiondata$position$x)

Nfilt <- 24
filt <- rep(1, Nfilt)/Nfilt

for (i.theta in 1:N.theta){
	##################################################
	## the true position
	tau <- theta.starts[i.theta]
	i.pos <- round(tau*10)
	ww.traj <- rbind(c(rep(1/33,33), rep(0, 67)), c(rep(0, 33), rep(1/34,34), rep(0, 33)), c(rep(0,67), rep(1/33, 33))) # weighted average for the represented position, variance; 100 Hz, 10 datapoints			

	##################################################
	## the posterior mean position
	mu.tau <- motiondata$post[[i.theta]]$estimpar$mu # for MAP, DDCal and PPC, the posterior in "motion-coordinates"
	L.theta <- ncol(mu.tau)
	
	T.theta <- motiondata$post[[i.theta]]$tstats[4] * 1000 # this is the real duration of the theta sequence in ms
	t.samples <- seq(1, T.theta, length=L.theta)
	t.samples.ms <- seq(1, T.theta)
	
	xy.ms <- rbind(approx(t.samples, mu.tau[1,], t.samples.ms)$y, approx(t.samples, mu.tau[2,], t.samples.ms)$y) # interpolating
	xy.ms <- cbind(matrix(rep(xy.ms[,1], Nfilt/2), 2), xy.ms, matrix(rep(xy.ms[,T.theta], Nfilt/2), 2)) # padding for filtering
	xy.ms <- rbind(filter(xy.ms[1,], filt), filter(xy.ms[2,], filt))[,(Nfilt/2):(T.theta+Nfilt/2-1)] # filtering

	##################################################
	## the posterior variance: txy_var[[i.theta]]
	## the represented position: 	txy_mu[[i.theta]]
	N_mu <- ncol(txy_mu[[i.theta]])
	L1 <- floor(N_mu / 3)
	L2 <- N_mu - 2*L1		
	ww.traj <- rbind(c(rep(1/L1,L1), rep(0, L1+L2)),  c(rep(0, L1), rep(1/L2,L2), rep(0, L1)),  c(rep(0,L1+L2), rep(1/L1, L1)))
	xy_theta_repr[,, i.theta] <- t(ww.traj %*% t(txy_mu[[i.theta]]))
	var_theta_post[,, i.theta] <- t(ww.traj %*% t(txy_var[[i.theta]]))
	mu_theta_post[,, i.theta] <- t(ww.traj %*% t(xy.ms))
	if (i.theta %% 1000 == 0) cat(i.theta, ' ')
}

	

################################################################
## CCV with different spike thresholds and mean vs. median
################################################################
nmax <- dim(decode.theta$postmean)[3]
ix <- 1:(nmax-1); iy <- 2:nmax

ths <- round(c(0, median(decode.theta$nsp)),1)
CCV_all <- list()
CCV_median <- list()

for (k in 2:1){
	sp_th <- ths[k]
	for (j in 1:3){
		if (j==1) i_supra_th <- (decode.theta$nsp[j,ix] > sp_th) & (decode.theta$nsp[j,iy] > sp_th)
		CCV_j <- colSums((decode.theta$postmean[,j,ix] - decode.theta$postmean[,j,iy])^2)
		CCV_jk <- CCV_j[i_supra_th]
		if (k==1) {
			CCV_all[[j]] <- CCV_jk
			TEE_all[[j]] <- TEE_all[[j]][c(i_supra_th, F)]	
		}
		if (k==2) {
			CCV_median[[j]] <- CCV_jk
			TEE_median[[j]] <- TEE_all[[j]][c(i_supra_th, F)]	
		}
	}
}

m_ccv <- array(NA, dim=c(2, 2, 3), dimnames=list(c('mean', 'median'), paste('threshold', ths), c('past', 'present', 'future')))
se_ccv <- array(NA, dim=c(2, 2, 3), dimnames=list(c('mean', 'median'), paste('threshold', ths), c('past', 'present', 'future')))




mbr <- ceiling(max(c(unlist(lapply(CCV_all, max)), unlist(lapply(TEE_all, max, na.rm=T)))))
brs <- seq(0, mbr, by=1/100)
mids <- brs[-1] - 1/200
CDFs <- array(NA, dim=c(4, 3, length(brs)-1), dimnames=list(c('TEE-all', 'TEE-hsc', 'CCV-all', 'CCV-hsc'), c('early', 'mid', 'late'), mids))
for (ii in 1:3){	
	CDFs[1, ii, ] <- cumsum(hist(TEE_all[[ii]], br=brs, plot=F)$density / 100)
	CDFs[2, ii, ] <- cumsum(hist(TEE_median[[ii]], br=brs, plot=F)$density / 100)
	CDFs[3, ii, ] <- cumsum(hist(CCV_all[[ii]], br=brs, plot=F)$density / 100)
	CDFs[4, ii, ] <- cumsum(hist(CCV_median[[ii]], br=brs, plot=F)$density / 100)
}
eml_cols <- c('#b9b9b9', '#ffbd40', '#ff7474')
	

pdf(file=paste('./SimFigs/Fig6/ccv_tee_', reg, '_', code, '_j', jitter, sampler, '.pdf', sep=''), 4, 4, useD=F)
par(mfcol=c(2, 2))
par(mar=c(2.5,3.8,2,1))

matplot(mids, t(CDFs[1,,]), lty=1, t='l', col=eml_cols, lwd=2, log='x', axes=F, ylim=c(0, 1), xlab='CCV and TEE', ylab='Cum. probability', main='all cycles', cex.lab=1.2)
axis(1, c(0.005, 0.05, 0.5, 5),  c('0.005', '0.05', '0.5', '5'), cex.axis=1.2)
axis(2, 0:2/2, las=2, cex.axis=1.2)
matplot(mids, t(CDFs[3,,]), lty=2, t='l', col=eml_cols, lwd=2, log='x', add=T)
legend('bottomright', c('early', 'mid', 'late'), lty=1, lwd=2, bty='n', col=eml_cols)
legend('topright', c('TEE', 'CCV'), lty=c(1,2), lwd=2, bty='n', col=eml_cols[1])

plot(1:3, lapply(CCV_all, mean, na.rm=T), ylim=c(0, 0.1), gap=0, aces=F, , cex=2, pch=21, bg=2, axes=F, xlab='', ylab='error (cm2)', t='o', xlim=c(0.5, 3.5), cex.lab=1.2)
plotCI(1:3, unlist(lapply(CCV_all, mean, na.rm=T)), unlist(lapply(CCV_all, se, na.rm=T)), gap=0, cex=2, pch=21, add=T)
points(1:3, lapply(TEE_all, mean, na.rm=T), cex=2, pch=23, bg=4, t='o')
plotCI(1:3, unlist(lapply(TEE_all, mean, na.rm=T)), unlist(lapply(TEE_all, se, na.rm=T)), pch=23, gap=0, add=T, cex=2)
axis(1, 1:3, c('early', 'mid', 'late'), tick=F, las=2, cex.axis=1.2)
legend('topleft', leg=c('ccv', 'tee'), pch=c(21,23), pt.bg=c(2,4), bty='n', pt.cex=1.2)
axis(2, 0:2/20, 0:2*500, las=2, cex.axis=1.2)

matplot(mids, t(CDFs[2,,]), lty=1, t='l', col=eml_cols, lwd=2, log='x', axes=F, ylim=c(0, 1), xlab='CCV and TEE', ylab='Cum. probability', main='high spike count', cex.lab=1.2)
axis(1, c(0.005, 0.05, 0.5, 5),  c('0.005', '0.05', '0.5', '5'), cex.axis=1.2)
axis(2, 0:2/2, las=2, cex.axis=1.2)
matplot(mids, t(CDFs[4,,]), lty=2, t='l', col=eml_cols, lwd=2, log='x', add=T)
legend('bottomright', c('early', 'mid', 'late'), lty=1, lwd=2, bty='n', col=eml_cols)

plot(1:3, lapply(CCV_median, mean, na.rm=T), ylim=c(0, 0.1), gap=0, aces=F, , cex=2, pch=21, bg=2, axes=F, xlab='', ylab='error (cm2)', t='o', xlim=c(0.5, 3.5), cex.lab=1.2)
plotCI(1:3, unlist(lapply(CCV_median, mean, na.rm=T)), unlist(lapply(CCV_median, se, na.rm=T)), gap=0, cex=2, pch=21, add=T)
points(1:3, lapply(TEE_median, mean, na.rm=T), cex=2, pch=23, bg=4, t='o')
plotCI(1:3, unlist(lapply(TEE_median, mean, na.rm=T)), unlist(lapply(TEE_median, se, na.rm=T)), pch=23, gap=0, add=T, cex=2)
axis(1, 1:3, c('early', 'mid', 'late'), tick=F, las=2, cex.axis=1.2)
legend('topleft', leg=c('ccv', 'tee'), pch=c(21,23), pt.bg=c(2,4), bty='n', pt.cex=1.2)
axis(2, 0:2/20, 0:2*500, las=2, cex.axis=1.2)
dev.off()

P_all <- t.test(CCV_all[[3]]-CCV_all[[1]] -(TEE_all[[3]] - TEE_all[[1]]))
P_median <- t.test(CCV_median[[3]]-CCV_median[[1]]-(TEE_median[[3]] - TEE_median[[1]]))


SI_all <- mean(CCV_all[[3]]-CCV_all[[1]] -(TEE_all[[3]] - TEE_all[[1]]), na.rm=T)
seSI_all <- sd(CCV_all[[3]]-CCV_all[[1]] -(TEE_all[[3]] - TEE_all[[1]]), na.rm=T) / sqrt(sum(!is.na(CCV_all[[3]] + CCV_all[[1]] + TEE_all[[3]] + TEE_all[[1]])))

SI_med <- mean(CCV_median[[3]]-CCV_median[[1]] -(TEE_median[[3]] - TEE_median[[1]]), na.rm=T)
seSI_med <- sd(CCV_median[[3]]-CCV_median[[1]] -(TEE_median[[3]] - TEE_median[[1]]), na.rm=T) / sqrt(sum(!is.na(CCV_median[[3]] + CCV_median[[1]] + TEE_median[[3]] + TEE_median[[1]])))

SI <- matrix(c(SI_all, seSI_all, P_all$p.value, SI_med, seSI_med, P_median$p.value), 3, dimnames=list(c('mean', 'SE', 'P'), c('all', 'high spike count')))

sampling_index <- list(SI=SI, CCV_all=CCV_all, CCV_median=CCV_median, TEE_all=TEE_all, TEE_median=TEE_median, P_all=P_all, P_median=P_median)
SIDfile <- paste('./SimFigs/Fig6/SamplingIndex_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, reg, '_', jitter, sampler, '.RData', sep='')
save(sampling_index, file= SIDfile)

