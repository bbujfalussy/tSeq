
require(circular)
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


decfile <- paste('R', i_rat, 'D', i_day, '_decoded_theta_seqs.RData', sep='')
tryCatch({load(file=decfile)}, error = function(e) print('File containing the decoded theta sequences not found. Run the MAP_decode.R first!'), warning=function(w) print('File containing the decoded theta sequences not found. Run the MAP_decode.R first!'))

load(file=decfile)
chunks <- dec_theta$chunks

estimpos_early <- dec_theta$est_early
estimpos_mid <- dec_theta$est_mid
estimpos_late <- dec_theta$est_late


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
all_trials[1,] <- c(theta_chunk[[1]]$pos[1,1]-0.1, all_trials[2,1])

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

endpoint_stats <- array(0, dim=c(1, 12), dimnames=list(NULL, c('distance', 'direction', 'end_x', 'end_y', 'nsp', 'type', 'chunk', 'tstart', 'xpos', 'ypos', 'xend', 'yend')))
#type : 1: home, 2: away, 0: neither

estimpos <- estimpos_late
theta_chunk <- chunks$future
				
for (i.theta in 1:length(estimpos)){
	epos <- estimpos[[i.theta]]
	t_chunk <- theta_chunk[[i.theta]]
	if (length(epos$dir) > 1){
		if (nrow(epos$dir) > 1){
			nsp_i <- rowSums(t_chunk$rast, na.rm=T)
			L <- length(nsp_i)
			i_type <- is_in_trials(t_chunk$pos[,1], home_trials) + 2 * is_in_trials(t_chunk$pos[,1], away_trials) # including home, away or all trials
 # including home, away or all trials
			for (i.cycle in 1:L){
				decoded_position <- matrix(epos$estim_pos[i.cycle,], 2, 1) / 100
				true_position <- matrix(t_chunk$pos[i.cycle,2:3], 2, 1) / 100
				motion_direction <- t_chunk$pos[i.cycle,4]
				tr_end <- rotateX(decoded_position, (-1) * motion_direction, origin=true_position, offset=F)
				dist_dir <- c(sqrt(sum(tr_end^2)), coord2rad(tr_end[1], tr_end[2]))

				endpoint_stats_cycle <- c(dist_dir, tr_end, nsp_i[i.cycle], i_type[i.cycle], i.theta, t_chunk$pos[i.cycle,1:3], decoded_position*100)
				endpoint_stats <- rbind(endpoint_stats, endpoint_stats_cycle)				
			}
			endpoint_stats_cycle <- c(rep(NA,4), 0, 0, i.theta, NA, NA, NA, NA, NA)
			endpoint_stats <- rbind(endpoint_stats, endpoint_stats_cycle)				
		}
	}
}

endpoint_stats <- endpoint_stats[-1,]
L <- nrow(endpoint_stats)

i.shuffle <- sample(1:L, L)
shuffle_endpoint_stats <- endpoint_stats[i.shuffle,1:2]

tra_d_dist <- diff(endpoint_stats[,1])
tra_d_dir <- apply(rbind(diff(endpoint_stats[,2]) %% (2*pi), rev(diff(rev(endpoint_stats[,2]))) %% (2*pi)), 2, min)
shuffle_d_dist <- diff(shuffle_endpoint_stats[,1])
shuffle_d_dir <- apply(rbind(diff(shuffle_endpoint_stats[,2]) %% (2*pi), rev(diff(rev(shuffle_endpoint_stats[,2]))) %% (2*pi)), 2, min)

ymax <- ceiling(10*max(c(hist(shuffle_d_dir, br=seq(0, pi, length=20), plot=F)$density, hist(tra_d_dir, br=seq(0, pi, length=20), plot=F)$density)))/10


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

LR <- (as.vector(sign(endpoint_stats[,2] - pi)) + 1) / 2
LR_shuffle <- (as.vector(sign(shuffle_endpoint_stats[,2] - pi)) + 1) / 2
CCP <- get_CCP(LR)
CCP_shuffle <- get_CCP(LR_shuffle)

CCP <- CCP[CCP>0]
CCP_shuffle <- CCP_shuffle[CCP_shuffle>0]




N_shuffle <- 10000
CCP4 <- matrix(NA, 2, N_shuffle + 1)
CCP4[1, N_shuffle + 1] <- sum(CCP > 3)
CCP4[2, N_shuffle + 1] <- sum(CCP > 3) / (sum(CCP==2) + sum(CCP==3))
for (i_shuffle in 1:N_shuffle){
	ii <- sample(1:L, L)
	shuffle_dir <- endpoint_stats[ii,2]
	LR_shuffle <- (as.vector(sign(shuffle_dir - pi)) + 1) / 2
	CCP_shuffle <- get_CCP(LR_shuffle)
	CCP4[1, i_shuffle] <- sum(CCP_shuffle > 3)
	CCP4[2, i_shuffle] <- sum(CCP_shuffle > 3) / (sum(CCP_shuffle==2) + sum(CCP_shuffle==3))
	if (i_shuffle %% 100 == 0) print(i_shuffle)
}

P_value <- sum(CCP4[1,1:N_shuffle] >= CCP4[1,1+N_shuffle]) / N_shuffle
P_value2 <- sum(CCP4[2,1:N_shuffle] >= CCP4[2,1+N_shuffle]) / N_shuffle

CCP <- CCP[CCP>0]
CCP_shuffle <- CCP_shuffle[CCP_shuffle>0]
max_br <- max(c(CCP, CCP_shuffle))
hCCP <- hist(CCP, br=seq(0, max_br)+0.5, plot=F)
hCCPs <- hist(CCP_shuffle, br=seq(0, max_br)+0.5, plot=F)

dir.create('./DataFigures/Fig7/', showW=F)

cycling_file <- paste('DataFigures/Fig7/R', i_rat, 'D', i_day, '_cycling_nolog.pdf', sep='')
pdf(file=cycling_file, 11, 3, useD=F)
ymax <- max(c(hCCP$counts, hCCPs$counts))
par(mfcol=c(1,3))

mp <- barplot(hCCPs$counts+1/10, ylim=c(10,ymax), col=grey(0.25), xlab='duration (theta cycles)', ylab='No. of cycles', main='cycling duration', log='', axes=F, space=0, xlim=c(0,30))
# mp <- barplot(hCCPs$counts+1/10, ylim=c(10,ymax), col=grey(0.25), xlab='duration (theta cycles)', ylab='No. of cycles', main='cycling duration', log='y', axes=F, space=0, xlim=c(0,30))
xx <- 1:13
px <- xx/(2^(xx+1))
px <- px / sum(px)
points(mp[1:13], px*L, col=2, pch=16, t='o')
axis(1)
axis(2, las=2)
barplot(hCCP$counts+1/10, col=viridis(11, alpha=0.75)[8], log='', axes=F, add=T, ylim=c(10,ymax), space=0, xlim=c(0,30))
# barplot(hCCP$counts+1/10, col=viridis(11, alpha=0.75)[8], log='y', axes=F, add=T, ylim=c(10,ymax), space=0, xlim=c(0,30))
legend('topright', leg=c('shuffle', 'data', 'expected'), bty='n', pch=c(22,22,21), pt.bg=c(grey(0.7), viridis(11, alpha=0.75)[8], 2), pt.cex=1.5)

ymin <- 100*floor(min(CCP4[1,]) / 100)
ymax <- 100*ceiling(max(CCP4[1,]) / 100)

maintext <- paste('theta cycles CCP > 3, P=', P_value, sep='')
hist(CCP4[1,1:N_shuffle], br=seq(ymin, ymax, by=10), xlab='No. of cycles', ylab='Shuffle count', main=maintext)
abline(v=CCP4[1,N_shuffle + 1], col=2)

ymin <- floor(min(CCP4[2,]) * 100) / 100
ymax <- ceiling(max(CCP4[2,]) * 100) / 100

maintext <- paste('theta cycles CCP > 3 / CCP(1,2), P=', P_value2, sep='')
hist(CCP4[2,1:N_shuffle], br=seq(ymin, ymax, by=0.01), xlab='No. of cycles', ylab='Shuffle count', main=maintext)
abline(v=CCP4[2,N_shuffle + 1], col=2)

dev.off()


CCPdata <- list(CCP=CCP, CCP_shuffle=CCP_shuffle, CCP4=CCP4)
CCP_file <- paste('DataFigures/Fig7/R', i_rat, 'D', i_day, '_CCP.RData', sep='')
save(CCPdata, file=CCP_file)

