########################################################
### Calculate CCP - the constant cycling periods for synthetic data
source('./ExpData/SpeedFunctions.R', chdir = TRUE)
require(circular)
require(viridis)
require(colormap)
library(scico)



############################################################################
### Helper functions
############################################################################

rotateX <- function(X, phi, origin=c(0,0), offset=F){
	# 2 x N matrix of position coordinates to rotate
	# phi: rotation angle in radians
	# origin: the center of the rotation
	# shift: if true, then the origin is added to the rotated coordinates
	phi <- as.vector(phi)
	RotM <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2, 2)
	X0 <- rbind(X[1,] - origin[1], X[2,] - origin[2])

	Xrot <- RotM %*% X0
	if (offset) {
		XX <- rbind(Xrot[1,] + origin[1], Xrot[2,] + origin[2])
	} else {
		XX <- Xrot
	}
	XX
}

sigm <- function(x, ampl=1, slope=1, threshold=0){
	fx <- ampl / (1 + exp(-slope*(x-threshold)))
	fx
}


get_CCP <- function(LR){ ## calculating length of continuous cycling periods
	## LR is a vector of -1 and 1 indicating whether the trajectory was on the left (-1) or right (1). 
	##      end of a theta oscillation is indicated with NA
	L <- length(LR)
	LR[which(is.na(LR))] <- 3 # contiguous candidate theta cycles (contig period) change is signaled by NA changed here to 3
	CCP <- rep(1, L) # length of continuous cycling period
	last_dir <- 3 # 3 means NA
	CCP_dur <- 1
	CCP_start <- 1

	for (i in 1:L){
		if (LR[i] < 3){
			if (last_dir > 1){ # start of contig period 
				CCP_start <- i
				CCP_dur <- 1
				last_dir <- LR[i]			
			} else {
				if (LR[i] == last_dir){ # end of CCP
					CCP[CCP_start:i] <- CCP_dur
					CCP_start <- i
					CCP_dur <- 1
					last_dir <- LR[i]
				} else { # CCP goes on
					CCP_dur <- CCP_dur + 1
					last_dir <- LR[i]
				}
			}
		} else { # NA direction: end of CCP
			CCP[CCP_start:i] <- CCP_dur
			CCP[i] <- 0
			CCP_start <- i+1
			CCP_dur <- 1		
			last_dir <- LR[i] # NA?
		}
	}
	CCP
}



############################################################################
### Loading the necessary data
############################################################################

Tmax <- 1800
N.cells <- 200
code <- 'sampling'
spPth <- 32
pars <- cbind(slope=c(-8, -5, 0, 8, 5, 5), th=c(pi/8, pi/4, pi/2, pi/8, pi/4, pi/2))
rownames(pars) <- list('strong +', 'weak +', 'indep', 'tiny -', 'weak -', 'strong -')
corr_types <- c('same', 'similar', 'uniform', 'anti-quarter', 'anti-half', 'strong-anti')

sim.irreg <- F
N.samples <- 100

motionfile <- paste('SynthData/motiondata_Tmax', Tmax, '.RData', sep='')
tryCatch({load(file=motionfile)}, error = function(e) print('error: motionfile not read'), warning=function(w) print('error: motionfile not found, run the sim_thetaseq_kalman.R script first'))

tt <- motiondata$T.axis
xpos <- motiondata$position$x[1,] * 100
ypos <- motiondata$position$x[2,] * 100
# motion_dir <- atan2(diff(ypos), diff(xpos))
# motion_dir <- c(motion_dir[1], motion_dir)
# pos <- cbind(tt, xpos, ypos, motion_dir) 
# txyvd <- posfilter(pos, sdfilt=0.25, graphics=T)
# stats <- plot_motion_stats(xyvd= txyvd[,2:5], dt_pos=dt.run)

tt_30Hz <- seq(1/30, Tmax, by=1/30)
xpos30 <- approx(tt, xpos, tt_30Hz, rule=2)$y
ypos30 <- approx(tt, ypos, tt_30Hz, rule=2)$y
motion_dir30 <- atan2(diff(ypos30), diff(xpos30))
motion_dir30 <- c(motion_dir30[1], motion_dir30)
pos30 <- cbind(tt_30Hz, xpos30, ypos30, motion_dir30) 
txyvd30 <- posfilter(pos30, sdfilt=0.25, graphics=F)
# stats <- plot_motion_stats(xyvd= txyvd30[,2:5], dt_pos=1/30, file='model_motion_stats_new.pdf')

taus <- motiondata$post[[1]]$tstat[1] # theta cycle times
for (tt in 2:length(motiondata$post)) taus <- c(taus, motiondata$post[[tt]]$tstat[1])


N.theta <- length(motiondata$post)
xpos.theta <- approx(txyvd30[,1], txyvd30[,2], taus, rule=2)$y
ypos.theta <- approx(txyvd30[,1], txyvd30[,3], taus, rule=2)$y
motion_dir.theta <- approx(txyvd30[,1], txyvd30[,5], taus, rule=2)$y


############################################################################
### Calculating the CCP
############################################################################

CCPdata <- matrix(NA, 6, 3, dimnames=list(NULL, c('shuffle_mean', 'shuffle_SD', 'data')))
CyclingIndexData <- matrix(NA, 6, 3, dimnames=list(NULL, c('shuffle_mean', 'shuffle_SD', 'data')))


## to calculate the CCP we can use two different reference points:
## 1. the true position and motion direction of the animal 
## 2. the estimated position and motion direction of the animal 
references <- c('inferred', 'true')

## similarly, we can use the encoded and the decoded trajectory endpoints here
endpoints <- c('encoded', 'decoded')
# note, that for real data we can only use the combination 'true' and 'decoded'

samcols <- scico(100, palette='broc')[c(30,40,50,60,70,80)]
samcols_alpha <- scico(100, palette='broc', alpha=0.4)[c(30,40,50,60,70,80)]
samborder <- c(rep(grey(0), 6))

dir.create('./SimFigs/', showW=F)
dir.create('./SimFigs/Fig7/', showW=F)

for (i_end in 1:2){
	endpoint <- endpoints[i_end]
	for (i_ref in 1:2){
		reference <- references[i_ref]
		filename <- paste('SynthData/SimFigs/Fig7/ConstantCyclingPeriods_', endpoint, '_', reference, '.pdf', sep='')
	
		pdf(file=filename, 15, 7, useD=F)
		par(mfcol=c(3,6))
		
		
		for (i_par in 1:6){
			decodefile <- paste('SynthData/decode_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, '_evenspikes_0_', corr_types[i_par], '.RData', sep='')
			tryCatch({load(file= decodefile)}, error = function(e) print('error: motionfile not read'), warning=function(w) print('error: motionfile not found, run the Decode_theta_sim.R script first'))
		
			datafile <- paste('SynthData/spikes_Tmax', Tmax, '_Ncells', N.cells, '_', code, '_spPth', spPth, '_', corr_types[i_par], '.RData', sep='')
			tryCatch({load(file= datafile)}, error = function(e) print('error: motionfile not read'), warning=function(w) print('error: motionfile not found, run the sim_thetaseq_kalman.R script first'))
			
			tra_dist_dir <- matrix(0, 2, N.theta)
			trajectory_ends_idio <- matrix(0, 2, N.theta)
			for (i.theta in 1:N.theta){	
				if (endpoint == 'encoded'){
					traj_end <- data$xy.trajectories[[i.theta]]
					kk <- 100
				}	else {
					traj_end <- decode.theta$postmean[,,i.theta]
					kk <- 3
				}
				# last point of the trajectory rotated around the TRUE position
				if (reference == 'true'){
					tr_end <- rotateX(traj_end, (-1) * motion_dir.theta[i.theta], origin=c(xpos.theta[i.theta]/100, ypos.theta[i.theta]/100), offset=F)[,kk]
				}
				
				# last point of the trajectory rotated around the INFERRED position
				if (reference == 'inferred'){
					dxy <- motiondata$post[[i.theta]]$smooth$mu[,11] - motiondata$post[[i.theta]]$smooth$mu[,10]
					inferred_direction <- coord2rad(dxy[1], dxy[2])
					tr_end <- rotateX(traj_end, (-1) * inferred_direction, origin=motiondata$post[[i.theta]]$smooth$mu[,11], offset=F)[,kk]
				}
				tra_dist_dir[,i.theta] <- c(sqrt(sum(tr_end^2)), coord2rad(tr_end[1], tr_end[2]))
				trajectory_ends_idio[,i.theta] <- tr_end
			}
			
			i.shuffle <- sample(1:N.theta, N.theta)
		
			LR <- (as.vector(sign(tra_dist_dir[2,] - pi)) + 1) / 2
			LR_shuffle <- (as.vector(sign(tra_dist_dir[2,i.shuffle] - pi)) + 1) / 2
			CCP <- get_CCP(LR)
			CCP_shuffle <- get_CCP(LR_shuffle)
		
			CCP <- CCP[CCP>0]
			CCP_shuffle <- CCP_shuffle[CCP_shuffle>0]
		
		
			N_shuffle <- 10000
			CCP4 <- matrix(NA, 2, N_shuffle + 1)
			CCP4[1, N_shuffle + 1] <- sum(CCP > 3) / length(CCP)
			CCP4[2, N_shuffle + 1] <- sum(CCP > 3) / (sum(CCP==2) + sum(CCP==3))
			for (i_shuffle in 1:N_shuffle){
				ii <- sample(1:N.theta, N.theta)
				shuffle_dir <- tra_dist_dir[2, ii]
				LR_shuffle <- (as.vector(sign(shuffle_dir - pi)) + 1) / 2
				CCP_shuffle <- get_CCP(LR_shuffle)
				CCP4[1, i_shuffle] <- sum(CCP_shuffle > 3) / length(CCP)
				CCP4[2, i_shuffle] <- sum(CCP_shuffle > 3) / (sum(CCP_shuffle==2) + sum(CCP_shuffle==3))
				if (i_shuffle %% 100 == 0) print(i_shuffle)
			}
		
			P_value <- sum(CCP4[1,1:N_shuffle] >= CCP4[1,1+N_shuffle]) / N_shuffle
			P_value2 <- sum(CCP4[2,1:N_shuffle] >= CCP4[2,1+N_shuffle]) / N_shuffle
		
			CCP <- CCP[CCP>0]
			CCP_shuffle <- CCP_shuffle[CCP_shuffle>0]
			# max_br <- max(c(CCP, CCP_shuffle))
			max_br <- 40
			hCCP <- hist(CCP, br=seq(0, max_br)+0.5, plot=F)
			hCCPs <- hist(CCP_shuffle, br=seq(0, max_br)+0.5, plot=F)
			
			# ymax <- max(c(hCCP$counts, hCCPs$counts))
			ymax <- 8000
			title_text <- rownames(pars)[i_par]	
			mp <- barplot(hCCPs$counts+1/10, space=0, log='', ylim=c(10,ymax), xlim=c(0, 20), col=grey(0.75), xlab='duration (theta cycles)', ylab='No. of cycles', main=title_text, axes=F)
			mp <- barplot(hCCP$counts+1/10, space=0, log='', col=samcols_alpha[i_par], add=T, ylim=c(10,ymax), xlim=c(0, 20), axes=F)
			axis(1, 0:8*5-0.5, 0:8*5)
			axis(2, c(10, 100, 1000, 5000), c(10, 100, 1000, 5000), las=2)
			xx <- 1:13
			px <- xx/(2^(xx+1))
			px <- px / sum(px)
			points(mp[1:13], px*N.theta, col=2, pch=16, t='o', cex=0.5)	
			legend('topright', leg=c('shuffle', 'data', 'expected'), bty='n', pch=c(22,22,21), pt.bg=c(grey(0.7), viridis(11, alpha=0.4)[8], 2), pt.cex=1.5)
			
			ymin <- floor(min(CCP4[1,]) * 100) / 100
			ymax <- ceiling(max(CCP4[1,]) * 100) / 100
			
			maintext <- paste('CCP > 3, P=', P_value, sep='')
			hist(CCP4[1,1:N_shuffle], br=seq(ymin, ymax, by=0.0025), xlab='prop. of cycles', ylab='Shuffle count', main=maintext)
			abline(v=CCP4[1,N_shuffle + 1], col=2)
		
			CCPdata[i_par,1] <- mean(CCP4[1,1:N_shuffle])
			CCPdata[i_par,2] <- sd(CCP4[1,1:N_shuffle])
			CCPdata[i_par,3] <- CCP4[1,1 + N_shuffle]	
			
			ymin <- floor(min(CCP4[2,]) * 100) / 100
			ymax <- ceiling(max(CCP4[2,]) * 100) / 100
			
			maintext <- paste('CCP > 3 / CCP(1,2), P=', P_value2, sep='')
			hist(CCP4[2,1:N_shuffle], br=seq(ymin, ymax, by=0.01), xlab='cycling index', ylab='Shuffle count', main=maintext)
			abline(v=CCP4[2,N_shuffle + 1], col=2)
			
			CyclingIndexData[i_par,1] <- mean(CCP4[2,1:N_shuffle])
			CyclingIndexData[i_par,2] <- sd(CCP4[2,1:N_shuffle])
			CyclingIndexData[i_par,3] <- CCP4[2,1 + N_shuffle]
			
		}	
		dev.off()	
	
		
		filename <- paste('SynthData/SimFigs/Fig7/ConstantCyclingPeriods_', endpoint, '_', reference, '_summary.pdf', sep='')
		pdf(file=filename, 5, 5, useD=F)
		par(mfcol=c(2,1))
		par(mar=c(4,4, 1,1))
		maintext <- 'CCP > 3'		
		if (endpoint == 'decoded'){
			ymin <- 0.15
			ymax <- 0.4
		} else {
			ymin <- 0
			ymax <- 0.8			
		}
		
		hist(CCP4[1,1:N_shuffle], br=seq(ymin, ymax, by=0.002), xlab='prop. of cycles', ylab='Shuffle count', main=maintext, axes=F)
		axis(1)
		axis(2, las=2)
		points(CCPdata[,3], rep(0, 6), pch=23, bg=samcols, lwd=1, cex=1.5)

		maintext <- 'CCP > 3 / CCP(2,3)'
		if (endpoint == 'decoded'){
			ymin <- 0.4
			ymax <- 1.4
		} else {
			ymin <- 0
			ymax <- 3.5	
		}
		
		hist(CCP4[2,1:N_shuffle], br=seq(ymin, ymax, by=0.008), xlab='Cycling index', ylab='Shuffle count', main=maintext, axes=F)
		axis(1)
		axis(2, las=2)
		points(CyclingIndexData[,3], rep(0, 6), pch=23, bg=samcols, lwd=1, cex=1.5)
		legend('topright', leg=rownames(pars), pch=23, pt.bg=samcols, bty='None')
		
		dev.off()	
		
		
	}	
	
}


