se <- function(x, ...) sqrt(var(x, ...)/sum(!is.na(x)))

is_in_trials <- function(times, start_end_matrix){
	# a function that checks for each value in the vector times 
	# whether it is within the intervals defined by the start_end_matrix
	# start_end_matrix is a matrix with 2 columns defining the start and the end of the target intervals
	#
	# return: a binary vector of length(times)
	is_in <- rep(F, length(times))
	k <- 1
	for (time in times){
		past_starts <- which(start_end_matrix[,1] < time)
		if (length(past_starts) > 0){
			i_row <- max(past_starts) # the last interval-start smaller than time	
			if (start_end_matrix[i_row,2] > time) is_in[k] <- T
		}
		k <- k + 1
	}
	is_in
}

decfile <- paste('R', i_rat, 'D', i_day, '_decoded_theta_seqs.RData', sep='')
if (decfile %in% list.files()) {
	print('loading decoded position')
	load(file=decfile)
	chunks <- dec_theta$chunks
	estimpos_early <- dec_theta$est_early
	estimpos_mid <- dec_theta$est_mid
	estimpos_late <- dec_theta$est_late
} else {
	print('decoded locations not found! Run MAP_decode.R first!')
}
	
dir.create('./DataFigures/Fig6', showW=F)

## calculating the number of spikes in theta cycles and the shift of the early, mid and late theta phase relative to the current time
quantile_curvature <- 0
savefile <- F
source('./decoding_curves_exp.R')

################################################################
## the difference between subsequent past, present and future representations
## standard MAP decoding

#################################################################################
## index of home and away trials
##################################################################################

if (i_day == 1) home_well <- 15 else home_well <- 29

wellinfo_file <- paste('Pfeiffer13/Rat', i_rat, 'Day', i_day, '_wellinfo.txt', sep='')
tryCatch({trls <- read.table(wellinfo_file, header = FALSE, sep = " ", dec = ".") }, error = function(e) print('File containing the reward locations not found. Please contact the authors to obtain it!'), warning=function(w) print('File containing the reward locations not found. Please contact the authors to obtain it!'))

trials <- rbind(trls, c(floor(tmax), NA, NA))
home_trials <- matrix(c(NA, NA), 1, dimnames=list(NULL, c('start', 'end')))
away_trials <- matrix(c(NA, NA), 1, dimnames=list(NULL, c('start', 'end')))
all_trials <- matrix(c(NA, NA), 1, dimnames=list(NULL, c('start', 'end')))
for (i_trial in 2:nrow(trials)){
	t_start <- trials[i_trial-1,1]
	t_end <- trials[i_trial,1]
	interval <- t_end - t_start
	if (trials[i_trial-1,2] != home_well) {
		away_trials <- rbind(away_trials, c(t_start, t_end))
		all_trials <- rbind(all_trials, c(t_start, t_end))
	} else {
		if (interval < 100){ # we exclude intervals where the animal did not find the goal within 100 s in home trials.
			home_trials <- rbind(home_trials, c(t_start, t_end))
			all_trials <- rbind(all_trials, c(t_start, t_end))
		}
	}
}
away_trials <- away_trials[-1,]
home_trials <- home_trials[-1,]
all_trials[1,] <- c(chunks$past[[1]]$pos[1,1]-0.1, all_trials[2,1])

if (include_trials == 'home'){
	trials_start_end <- home_trials
} else if (include_trials == 'away'){
	trials_start_end <- away_trials
} else {
	trials_start_end <- all_trials	
}
################################################################################
## TEE and CCV should be calculated from the same theta cycles
## we construct two lists with CCV and TEE in early, mid and late theta spikes

thresholds <- c(0, median(nsp))


for (i_th in 1:2){
	sp_th <- thresholds[i_th]
	i_supra_list  <- list()

	ccv_eml <- list() # eml: Early, Mid and Late
	tee_eml <- list()
	
	for (i.eml in 1:3){
		if (i.eml == 1) {
			estimpos <- estimpos_early
			theta_chunk <- chunks$past
		} else if (i.eml == 2){
			estimpos <- estimpos_mid
			theta_chunk <- chunks$present
		} else {
			estimpos <- estimpos_late
			theta_chunk <- chunks$future
		}		
		
		ccv_theta <- 0
		tee_theta <- 0
				
		for (i.theta in 1:length(estimpos)){
			epos <- estimpos[[i.theta]]
			t_chunk <- theta_chunk[[i.theta]]
			if (length(epos$dir) > 1){
				if (nrow(epos$dir) > 1){
					# theta cycles that pass the threshold - note, that both cycle i and cycle i+1 have to pass the threshold!
					if (i.eml ==1){
						nsp_i <- rowSums(t_chunk$rast, na.rm=T)
						L <- length(nsp_i)
						i_supra_th <- ((nsp_i[1:(L-1)] >= sp_th) & (nsp_i[2:L] >= sp_th)) # high number of spikes...

						i_type <- is_in_trials(t_chunk$pos[,1], trials_start_end)[1:L-1] # including home, away or all trials
						i_supra_type <- i_supra_th & i_type

						i_supra_list[[i.theta]] <- i_supra_type
					} else { # we only check in the early phase
						i_supra_type <- i_supra_list[[i.theta]]
					}
					ccv <- diff(epos$estim_pos[,1]/100)^2 + diff(epos$estim_pos[,2]/100)^2
					
					# reference position:
					xx_ref <- approx(txyvd[,1], txyvd[,2], xout=t_chunk$pos[,1] + delay_eml[i.eml])$y
					yy_ref <- approx(txyvd[,1], txyvd[,3], xout=t_chunk$pos[,1] + delay_eml[i.eml])$y
					tee <- rowSums((cbind(xx_ref, yy_ref)/100 - epos$estim_pos/100)^2)
					
					ccv_theta <- c(ccv_theta , ccv[i_supra_type])
					tee_theta <- c(tee_theta, tee[c(i_supra_type, F)])
					
				}
			}
		}
		ccv_eml[[i.eml]]  <- ccv_theta[-1]
		tee_eml[[i.eml]]  <- tee_theta[-1]
	}
	if (i_th == 1){
		CCV_all <- rbind(ccv_eml[[1]], ccv_eml[[2]], ccv_eml[[3]])
		TEE_all <- 	rbind(tee_eml[[1]], tee_eml[[2]], tee_eml[[3]])
	} else {
		CCV_median <- 	rbind(ccv_eml[[1]], ccv_eml[[2]], ccv_eml[[3]])
		TEE_median <- 	rbind(tee_eml[[1]], tee_eml[[2]], tee_eml[[3]])		
	}
}

cols <- c('#b9b9b9', '#ffbd40', '#ff7474')

mbr <- ceiling(max(c(max(unlist(CCV_all), na.rm=T), max(unlist(TEE_all), na.rm=T))))
brs <- seq(0, mbr, by=1/100)
mids <- brs[-1] - 1/200
CDFs <- array(NA, dim=c(2, 3, length(brs)-1), dimnames=list(c('TEE', 'CCV'), c('early', 'mid', 'late'), mids))
CDFs_median <- array(NA, dim=c(2, 3, length(brs)-1), dimnames=list(c('TEE', 'CCV'), c('early', 'mid', 'late'), mids))
for (ii in 1:3){
	CDFs_median[1, ii, ] <- cumsum(hist(TEE_median[ii,], br=brs, plot=F)$density / 100)
	CDFs_median[2, ii, ] <- cumsum(hist(CCV_median[ii,], br=brs, plot=F)$density / 100)		
	CDFs[1, ii, ] <- cumsum(hist(TEE_all[ii,], br=brs, plot=F)$density / 100)
	CDFs[2, ii, ] <- cumsum(hist(CCV_all[ii,], br=brs, plot=F)$density / 100)
}


mCCV_all <- apply(CCV_all, 1, mean)
mCCV_high <- apply(CCV_median, 1, mean)
seCCV_all <- apply(CCV_all, 1, se)
seCCV_high <- apply(CCV_median, 1, se)

mTEE_all <- apply(TEE_all, 1, mean)
mTEE_high <- apply(TEE_median, 1, mean)
seTEE_all <- apply(TEE_all, 1, se)
seTEE_high <- apply(TEE_median, 1, se)

### bootstrapping for confidence intervals on median
L <- length(CCV_all[1,])
Lmed <- length(CCV_median[3,])
N <- round(L/2)
Nmed <- round(Lmed/2)
K <- 1000
TEE_bootstrap <- matrix(NA, 3, K)
TEEhigh_bootstrap <- matrix(NA, 3, K)
CCV_bootstrap <- matrix(NA, 3, K)
CCVhigh_bootstrap <- matrix(NA, 3, K)
SI_bootstrap <- rep(NA, K)
SIhigh_bootstrap <- rep(NA, K)

for (k in 1:K){
	ii <- sample(1:L, N)
	TEE_bootstrap[,k] <- apply(TEE_all[,ii],1, median, na.rm=T)
	CCV_bootstrap[,k] <- apply(CCV_all[,ii],1, median, na.rm=T)
	SI_bootstrap[k] <- median(CCV_all[3,ii]-CCV_all[1,ii] - TEE_all[3,ii] + TEE_all[1,ii])
	
	ii <- sample(1:Lmed, Nmed)
	TEEhigh_bootstrap[,k] <- apply(TEE_median[,ii],1, median, na.rm=T)
	CCVhigh_bootstrap[,k] <- apply(CCV_median[,ii],1, median, na.rm=T)
	SIhigh_bootstrap[k] <- median(CCV_median[3,ii]-CCV_median[1,ii] - TEE_median[3,ii] + TEE_median[1,ii])
}

medCCV_all <- apply(CCV_all, 1, median)
medCCV_high <- apply(CCV_median, 1, median)
medTEE_all <- apply(TEE_all, 1, median)
medTEE_high <- apply(TEE_median, 1, median)

CI_CCV_all <- rbind(apply(CCV_bootstrap, 1, quantile, 0.05), apply(CCV_bootstrap, 1, quantile, 0.95))
CI_CCV_high <- rbind(apply(CCVhigh_bootstrap, 1, quantile, 0.05), apply(CCVhigh_bootstrap, 1, quantile, 0.95))

CI_TEE_all <- rbind(apply(TEE_bootstrap, 1, quantile, 0.05), apply(TEE_bootstrap, 1, quantile, 0.95))
CI_TEE_high <- rbind(apply(TEEhigh_bootstrap, 1, quantile, 0.05), apply(TEEhigh_bootstrap, 1, quantile, 0.95))


discfile <- paste('DataFigures/Fig6/R', i_rat, 'D', i_day, '_SI_index_', include_trials, '.pdf', sep='')
pdf(file= discfile, 6, 6, useD=F)

## ALL SPIKES
par(mfrow=c(3,3))
par(mar=c(3,5,3,1))
matplot(mids, t(CDFs[1,,]), lty=1, t='l', col=cols, lwd=2, log='x', axes=F, ylim=c(0, 1), xlab='error (m2)', ylab='cumulativ probability', main='all cycles', cex.lab=1.5)
axis(1, c(0.005, 0.05, 0.5, 5), c('0.005', '0.05', '0.5', '5'), cex.axis=1.5)
axis(2, las=2, c(0, 0.5, 1), cex.axis=1.5)
matplot(mids, t(CDFs[2,,]), lty=2, t='l', col=cols, lwd=2, log='x', add=T)
legend('bottomright', c('early', 'mid', 'late'), lty=1, lwd=2, bty='n', col=cols)
legend('topright', c('TEE', 'CCV'), lty=c(1, 2), lwd=2, bty='n', col=cols[1])

plot(1:3, mCCV_all, pch=21, bg=viridis(5)[3], ylab='mean error (cm2)', axes=F, main='all', xlab='', t='o', cex=2, xlim=c(0.5, 3.5), ylim=c(0, 0.2), cex.lab=1.5)
plotCI(1:3, mCCV_all, seCCV_all, gap=0, add=T, cex=2)
points(1:3, mTEE_all, pch=23, bg=grey(0.75), cex=2, t='o')
plotCI(1:3, mTEE_all, seTEE_all, gap=0, add=T, cex=2, pch=5)
axis(2, 0:2/10, 0:2*1000, las=2, cex.axis=1.5)
axis(1, c(1,2,3), c('early', 'mid', 'late'), tick=F, cex.axis=1.5, las=2)

plot(1:3, medCCV_all, pch=21, bg=viridis(5)[3], ylab='median error (cm2)', axes=F, main='all', xlab='', t='o', cex=2, xlim=c(0.5, 3.5), ylim=c(0, 0.2), cex.lab=1.5)
plotCI(1:3, medCCV_all, li=CI_CCV_all[1,], ui=CI_CCV_all[2,], gap=0, add=T, cex=2)
points(1:3, medTEE_all, pch=23, bg=grey(0.75), cex=2, t='o')
plotCI(1:3, medTEE_all, li=CI_TEE_all[1,], ui=CI_TEE_all[2,], gap=0, add=T, cex=2, pch=5)
axis(2, 0:2/10, 0:2*1000, las=2, cex.axis=1.5)
axis(1, c(1,2,3), c('early', 'mid', 'late'), tick=F, cex.axis=1.5, las=2)

## ABOVE MEDIAN - High spike couint cycles
matplot(mids, t(CDFs_median[1,,]), lty=1, t='l', col=cols, lwd=2, log='x', axes=F, ylim=c(0, 1), xlab='error (m2)', ylab='cumulartive probability', main='high spike count cycles', cex.lab=1.5)
axis(1, c(0.005, 0.05, 0.5, 5), c('0.005', '0.05', '0.5', '5'), cex.axis=1.5)
axis(2, las=2, c(0, 0.5, 1), cex.axis=1.5)
matplot(mids, t(CDFs_median[2,,]), lty=2, t='l', col=cols, lwd=2, log='x', add=T)
legend('bottomright', c('early', 'mid', 'late'), lty=1, lwd=2, bty='n', col=cols)

plot(1:3, mCCV_high, pch=21, bg=viridis(5)[3], ylab='mean error (cm2)', axes=F, main='high', xlab='', t='o', cex=2, xlim=c(0.5, 3.5), ylim=c(0, 0.1), cex.lab=1.5)
plotCI(1:3, mCCV_high, seCCV_high, gap=0, add=T, cex=2)
points(1:3, mTEE_high, pch=23, bg=grey(0.75), cex=2, t='o')
plotCI(1:3, mTEE_high, seTEE_high, gap=0, add=T, cex=2, pch=5)
axis(2, 0:2/20, 0:2*500, las=2, cex.axis=1.5)
axis(1, c(1,2,3), c('early', 'mid', 'late'), tick=F, cex.axis=1.5, las=2)

plot(1:3, medCCV_high, pch=21, bg=viridis(5)[3], ylab='median error (cm2)', axes=F, main='high', xlab='', t='o', cex=2, xlim=c(0.5, 3.5), ylim=c(0, 0.1), cex.lab=1.5)
plotCI(1:3, medCCV_high, li=CI_CCV_high[1,], ui=CI_CCV_high[2,], gap=0, add=T, cex=2)
points(1:3, medTEE_high, pch=23, bg=grey(0.75), cex=2, t='o')
plotCI(1:3, medTEE_high, li=CI_TEE_high[1,], ui=CI_TEE_high[2,], gap=0, add=T, cex=2, pch=5)
axis(2, 0:2/20, 0:2*500, las=2, cex.axis=1.5)
axis(1, c(1,2,3), c('early', 'mid', 'late'), tick=F, cex.axis=1.5, las=2)



SI_all <- CCV_all[3,]-CCV_all[1,] - TEE_all[3,] + TEE_all[1,]
SI_high <- CCV_median[3,] - CCV_median[1,] - TEE_median[3,] + TEE_median[1,]

SI_mat <- c(mean(SI_all), median(SI_all), mean(SI_high), median(SI_high))
uiSI_mat <- c(mean(SI_all)+se(SI_all), quantile(SI_bootstrap, 0.95), mean(SI_high) + se(SI_high), quantile(SIhigh_bootstrap, 0.95))
liSI_mat <- c(mean(SI_all)-se(SI_all), quantile(SI_bootstrap, 0.05),  mean(SI_high) - se(SI_high), quantile(SIhigh_bootstrap, 0.05))

barplot2(SI_mat, plot.ci = T, ci.l = liSI_mat, ci.u = uiSI_mat, col=c(grey(0.75), grey(0.85)), axes=F, ylim=c(0, 0.06), names=c('all-mean', 'all-median', 'hsp-mean', 'hsp-median'), las=2)
axis(2, 0:3/50, 0:3/50*10000, las=2)

P_all <- t.test(SI_all)
P_high <- t.test(SI_high)

dev.off()

sampling_index <- list(SI_all= SI_all, SI_high= SI_high, SI_mat=SI_mat, uiSI_mat= uiSI_mat, liSI_mat= liSI_mat)
sampling_file <- paste('DataFigures/Fig6/R', i_rat, 'D', i_day, '_sampling_index_', include_trials, '.RData', sep='')
save(sampling_index, file= sampling_file)

