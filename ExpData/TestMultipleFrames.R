

library(R.matlab)
library(colormap)
library(lattice)
library(scico)
library(ggplot2)


source('SpeedFunctions.R')
source('PlaceCellFunctions.R')
# source('FunGaussSpikes.R')
source('../SynthData/Utils.R')

dir.create('./DataFigures/Fig6_FS2', showW=F)

######################################################


sampling.rate <- 1250 # Hz
tryCatch({data <- readMat("Pfeiffer13/spikes.mat")}, error = function(e) print('File containing the spikes not found. Please contact the authors to obtain it!'), warning=function(w) print('File containing the spikes not found. Please contact the authors to obtain it!'))
data <- data[[1]]
names(data) <- list('rat1', 'rat2', 'rat3', 'rat4')
summary(data)
for (i in 1:4){
	names(data[[i]]) <- list('day1', 'day2')
	for (j in 1:2) names(data[[i]][[j]]) <- list('spikes', 'position', 'I_Cells', 'E_Cells')	
} 

#####################################
## process the position data

i_rat <- 1
i_day <- 1


if ('cells_remap.Rdata' %in% list.files('DataFigures/Fig6_FS2/')) {
	load(file='DataFigures/Fig6_FS2/cells_remap.Rdata')
	head(cells_remap)
} else {
	print ('calculating remapping statistics for each session. This takes about 1h.')

	cells_remap <- data.frame(matrix(0, 1, 5, dimnames=list(NULL, c('rat-day', 'z-difference', 'cross_norm-difference', 'significance', 'rate'))))
	# normalised magnitude: 0: same as between random splits; 1: same as between independent cells; 
	# significance: 0: not, 1: significant, -1: not reliable (independent splits are as different as different cells)

	for (i_rat in 1:4){
		for (i_day in 1:2){
			print(paste('rat: ', i_rat, ', day:', i_day, sep=''))
			if (i_day == 1) home_well <- 15 else home_well <- 29
			ddata <- data[[i_rat]][[i_day]]
			ncells <- max(ddata$spikes[,2])
			
			## filter position to remove noise and make it smoother
			txyvd <- posfilter(ddata$position, sdfilt=0.25, graphics=F)
			
			
			## detect periods of movement - at least 1s long move with at least 0.2 s stops between
			dt_pos <- median(diff(txyvd[,1]))
			vmin <- 5 # cm/s
			T_run_min <- 1 # s
			T_stop_min <- 0.2 # s
			ii_run <- txyvd[,4] > vmin
			ii_run[is.na(ii_run)] <- F
			i_run <- get_long_periods(ii_run, T_run_min / dt_pos, T_stop_min / dt_pos)
			
			## change the speed vector according to the run periods
			vv <- txyvd[,4]
			vv[i_run & (vv<vmin)] <- vmin + 1/10
			vv[!(i_run) & (vv>vmin)] <- vmin - 1/10
			
			# xyvd <- txyvd[,2:5]
			# xyvd[!i_run] <- NA
					
			#################################################################################
			## control ratemaps
			##################################################################################
			xy_brs <- seq(0, 200, by=5) # seq(-2.5, 202.5, by=5)
			posmap <- estimate_posmap(txyvd[,1:3], i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=T)
	
			spt <- ddata$spikes
			colnames(spt) <- c('time', 'cell')
			
			tmin <- max(min(txyvd[,1]), min(spt[,1]))
			tmax <- min(max(txyvd[,1]), max(spt[,1]))
			rmin <- 1/(tmax - tmin)
			
			ncells <- max(ddata$spikes[,2])
			ratemaps <- array(NA, dim=c(ncells, 40, 40))
			nspikes <- rep(NA, ncells)
			
			for (ii in 1:ncells){
				spi <- spt[spt[,2] == ii,1]
				ratemap <- estimate_ratemap(spi, txyvd[,1:3], i_data=i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap, graphics=0)
				ratemaps[ii,,] <- ratemap
				nspikes[ii] <- attr(ratemap, 'nsp')
			}
	
			attr(ratemaps, 'xcenters') <- attr(ratemap, 'xcenters') 
			attr(ratemaps, 'ycenters') <- attr(ratemap, 'ycenters') 
			
			mean_rate <- apply(ratemaps, 1, mean)
			sd_rate <- apply(ratemaps, 1, sd)
	
			#################################################################################
			## ratemaps from home and away trials
			##################################################################################
	
			wellinfo_file <- paste('Pfeiffer13/Rat', i_rat, 'Day', i_day, '_wellinfo.txt', sep='')
			tryCatch({trls <- read.table(wellinfo_file, header = FALSE, sep = " ", dec = ".") }, error = function(e) print('File containing the reward locations not found. Please contact the authors to obtain it!'), warning=function(w) print('File containing the reward locations not found. Please contact the authors to obtain it!'))
	
			trials <- rbind(trls, c(floor(tmax), NA, NA))
			ii_home <- rep(F, length(ii_run))
			ii_away <- rep(F, length(ii_run))
			home_intervals <- c()
			away_intervals <- c()
			for (i_trial in 2:nrow(trials)){
				i_start <- which.min((trials[i_trial-1,1] - txyvd[,1])^2)
				i_end <- which.min((trials[i_trial,1] - txyvd[,1])^2)
				interval <- trials[i_trial,1] - trials[i_trial-1,1]
				if (trials[i_trial-1,2] != home_well) {
					ii_away[i_start:i_end] <- T
					away_intervals <- c(away_intervals, interval)	
				} else {
					if (interval < 100){ # we exclude intervals where the animal did not find the goal within 100 s in home trials.
						ii_home[i_start:i_end] <- T 
						home_intervals <- c(home_intervals, interval)	
					}
				}
			}
	
			txyvdh <- cbind(txyvd, rep(0, nrow(txyvd)))
			colnames(txyvdh)[6] <- 'frame'
			ratdata <- data.frame(txyvdh)
			ratdata[,6] <- 'rest'
			ratdata[ii_away&i_run,6] <- 'away'
			ratdata[ii_home&i_run,6] <- 'home'
			ratdata$frame <- as.factor(ratdata$frame)
			ratdata <- ratdata[!is.na(ratdata$speed2),]
			maintext=paste('log P = ', round(log10(t.test(txyvd[ii_home&i_run, 4], txyvd[ii_away&i_run, 4])$p.value)))
	
			SpeedDiff_plot_file = paste('DataFigures/Fig6_FS2/ViolSpeed_rat', i_rat, '_day', i_day, '.pdf', sep='')
			pdf(file=SpeedDiff_plot_file, 5, 3, useD=F)		
			p <-ggplot(ratdata, aes(x=frame, y=speed2, col=frame), main=maintext) + geom_violin()		
			p + stat_summary(fun.data='mean_sdl', na.rm=T, fun.args = list(mult = 1), geom='pointrange', color=1) + ggtitle(maintext) + theme_minimal() + ylim(-10, 80)
			dev.off()
	
			posmap_home <- estimate_posmap(txyvd[,1:3], i_run&ii_home, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F)
			posmap_away <- estimate_posmap(txyvd[,1:3], i_run&ii_away, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F)
	
			## we normalise ratemaps by the global mean and SD...
			ratemaps_home <- array(NA, dim=c(ncells, 40, 40))
			ratemaps_away <- array(NA, dim=c(ncells, 40, 40))
	
			for (ii in 1:ncells){
				spi <- spt[spt[,2] == ii,1]
				ratemap <- estimate_ratemap(spi, txyvd[,1:3], i_data=i_run&ii_home, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap_home, graphics=0)
				ratemaps_home[ii,,] <- (ratemap - mean_rate[ii]) / sd_rate[ii]
	
				ratemap <- estimate_ratemap(spi, txyvd[,1:3], i_data=i_run&ii_away, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap_away, graphics=0)
				ratemaps_away[ii,,] <- (ratemap - mean_rate[ii]) / sd_rate[ii]
			}
	
	
			#########################################################
			# # random splitting the data to home and away trials
			#########################################################
	
			FramesTest_file <- paste('Rat', i_rat, 'Day', i_day, '_framesTest.RData', sep='')
			if (FramesTest_file %in% list.files('DataFigures/Fig6_FS2/')) {
				FramesTest_file <- paste('DataFigures/Fig6_FS2/', FramesTest_file, sep='')
				load(FramesTest_file)
			} else {
				K <- 100 # number of random splits
				maxN <- min(100, ncells-1)
		
				diff_away_home <- array(NA, dim=c(ncells, K+101))
				diff_away_home[,K+1] <- apply((ratemaps_home - ratemaps_away)^2, 1, mean)
				
				maxN <- min(100, ncells-1)
				for (ii in 1:ncells){			
					testset <- sample(setdiff(1:ncells, ii), maxN, replace=F)
					for (k in 1: maxN){
						jj <- testset[k]
						diff_away_home[ii,K+1+k] <- mean((ratemaps_home[ii,,] - ratemaps_away[jj,,])^2)				
					}
				}
		
		
				for (k in 1:K){
					print(paste('calculating random split #', k, ' ...', sep=''))
		
					### we sample randomly from the away or the home intervals
					# randtrials <- trials
					trial_start <- trials[1,1]
					i_trial <- 1
					randtrials <- array(c(trial_start, home_well, i_trial), dim=c(1,3))
					while (trial_start < tmax){
						i_trial <- i_trial + 1
						if (randtrials[i_trial-1,2] == home_well){ # previous goal was the home well
							well <- home_well + 1
							trial_start <- trial_start + sample(home_intervals, 1)
						} else {  # goal was the previous home well
							well <- home_well
							trial_start <- trial_start + sample(away_intervals, 1)					
						}
						randtrials <- rbind(randtrials, c(trial_start, well, i_trial))
					}
					randtrials <- randtrials[1:(i_trial-1),]
		
					# randtrials[,1] = sort(runif(nrow(trials), min(trials[,1]), max(trials[,1])))
					ii_home <- rep(F, length(ii_run))
					ii_away <- rep(F, length(ii_run))
					for (i_trial in 2:nrow(randtrials)){
						i_start <- which.min((randtrials[i_trial-1,1] - txyvd[,1])^2)
						i_end <- which.min((randtrials[i_trial,1] - txyvd[,1])^2)
						if (randtrials[i_trial-1,2] == home_well) ii_home[i_start:i_end] <- T else ii_away[i_start:i_end] <- T
					}
			
					posmap_home_rand <- estimate_posmap(txyvd[,1:3], i_run&ii_home, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F)
					posmap_away_rand <- estimate_posmap(txyvd[,1:3], i_run&ii_away, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F)
			
					## we normalise ratemaps by the global mean and SD...
					ratemaps_home_rand <- array(NA, dim=c(ncells, 40, 40))
					ratemaps_away_rand <- array(NA, dim=c(ncells, 40, 40))
			
					for (ii in 1:ncells){
						spi <- spt[spt[,2] == ii,1]
						ratemap <- estimate_ratemap(spi, txyvd[,1:3], i_data=i_run&ii_home, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap_home_rand, graphics=0)
						ratemaps_home_rand[ii,,] <- (ratemap - mean_rate[ii]) / sd_rate[ii]
			
						ratemap <- estimate_ratemap(spi, txyvd[,1:3], i_data=i_run&ii_away, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap_away_rand, graphics=0)
						ratemaps_away_rand[ii,,] <- (ratemap - mean_rate[ii]) / sd_rate[ii]
			
					}
					diff_away_home[,k] <- apply((ratemaps_home_rand - ratemaps_away_rand)^2, 1, mean)
				}
		
				FramesTest_file <- paste('DataFigures/Fig6_FS2/Rat', i_rat, 'Day', i_day, '_framesTest.RData', sep='')
				save(diff_away_home, file=FramesTest_file)
			}
			
			# 1. z-scoring the data with the distribution of the PFdiff across random splits
			PFdiff <- diff_away_home
			for (ii in 1:ncells){
				if (mean_rate[ii] > 0.1){
					m_diff <- mean(diff_away_home[ii,1:100])
					sd_diff <- sd(diff_away_home[ii,1:100])
					# z-scoring
					PFdiff[ii, ] <- (diff_away_home[ii,] - m_diff) / sd_diff
				} else {
					PFdiff[ii, ] <- 0
				}
			}
			
			# 2. We should focus on reliable cells: where the PFdiff across cells is larger than the PFdiff across random splits of the data. 
			## calculate the mean PFdiff across cells
			PFdiff_across_cells <- apply(PFdiff[,102:201], 1, mean, na.rm=T)
			ii_reliable <- PFdiff_across_cells > 2
			sum(ii_reliable)
			
			
			## Find the cells with significant PFdiff in home versus away
			PFdiff_significant <- rep(F, ncells)
			for (i_cell in 1:ncells){
				if ((PFdiff[i_cell,101] > sort(PFdiff[i_cell,1:100])[95]) & ii_reliable[i_cell]) PFdiff_significant[i_cell] <- T
			}
			
			cells_remap_rat <- data.frame(matrix(0, ncells, 5, dimnames=list(NULL, c('rat-day', 'z-difference', 'cross_norm-difference', 'significance', 'rate'))))
			cells_remap_rat[,1] <- paste('rat', i_rat, 'day', i_day, sep='')
			cells_remap_rat[!ii_reliable,4] <- -1
			cells_remap_rat[PFdiff_significant,4] <- 1
			cells_remap_rat[,5] <- mean_rate
			
			for (i_cell in 1:ncells){
				if (mean_rate[ii] > 0.1){
					cells_remap_rat[i_cell,2] <- PFdiff[i_cell,101]
					mean_cross_diff <- mean(PFdiff[i_cell,1:100+101], na.rm=T)
					cells_remap_rat[i_cell,3] <- PFdiff[i_cell,101] / mean_cross_diff
				} else {
					cells_remap_rat[i_cell,2] <- 0
					cells_remap_rat[i_cell,3] <- 0
				}
			}
			head(cells_remap_rat)
			cells_remap <- rbind(cells_remap, cells_remap_rat)
	
					
			PFdiff_plot_file = paste('DataFigures/Fig6_FS2/PFields_rat', i_rat, '_day', i_day, '.png', sep='')
			Nmax <- 12
			png(file=PFdiff_plot_file, Nmax*100, 300)
			
			par(mfcol=c(3,Nmax))
			par(mar=c(1,1,1,1)/5)
			for (ii in which(PFdiff_significant&ii_reliable)[1:Nmax+23]){
				rr <- ratemaps[ii,,]
				rr_h <- ratemaps_home[ii,,] * sd_rate[ii] + mean_rate[ii]			
				rr_a <- ratemaps_away[ii,,] * sd_rate[ii] + mean_rate[ii]
				rr[rr < 0] <- 0
				rr_a[rr_a < 0] <- 0
				rr_h[rr_h < 0] <- 0
				maxrate <-max(c(max(rr), max(rr_h), max(rr_a)))
				brs = seq(0, ceiling(maxrate), length=101)
				maintext <- round(cells_remap_rat$cross_norm.difference[ii],2)
				image(rr, col=viridis(100, option='B'), axes=F, br=brs)
				text(0.5, 1, round(max(rr)), col= viridis(100, option='B')[99], pos=1, cex=2)
				text(0.5, 0, maintext, col= viridis(100, option='B')[99], pos=3, cex=2)
				image(rr_h, col=viridis(100, option='B'), axes=F, br=brs)
				text(0.5, 1, round(max(rr_h)), col= viridis(100, option='B')[99], pos=1, cex=2)
				image(rr_a, col=viridis(100, option='B'), axes=F, br=brs)
				text(0.5, 1, round(max(rr_a)), col= viridis(100, option='B')[99], pos=1, cex=2)
				
			}
			dev.off()		
								
					
			mean_rate <- apply(ratemaps, 1, mean)
			mean_rate_h <- mean_rate
			mean_rate_a <- mean_rate
	
			max_rate <- apply(ratemaps, 1, max)
			max_rate_h <- mean_rate
			max_rate_a <- mean_rate
			
			for (ii in 1:ncells){
				rr_h <- ratemaps_home[ii,,] * sd_rate[ii] + mean_rate[ii]
				mean_rate_h[ii] <- mean(rr_h)
				max_rate_h[ii] <- max(rr_h)
				rr_a <- ratemaps_away[ii,,] * sd_rate[ii] + mean_rate[ii]
				max_rate_a[ii] <- max(rr_a)
			}
				
		}
	}	
	
	cells_remap <- cells_remap[-1,]
	cells_remap$rat.day <- as.factor(cells_remap$rat.day)
	save(cells_remap, file='DataFigures/Fig6_FS2/cells_remap.Rdata')
}

#######################################################################
## percent of significant cells (of the reliable spatially tuned cells)
percent_remap <- rep(NA, 8)
for (i_data in 1:8){
	i_session <- cells_remap$rat.day == levels(cells_remap$rat.day)[i_data]
	N_reliable <- sum(cells_remap$significance[i_session] > -1)
	N_remap <- sum(cells_remap$significance[i_session] > 0)
	percent_remap[i_data] <- N_remap/N_reliable
}

ratcols <- c("#FF2151", "#FF8FA8", "#DF00B3", "#FF7FF5", "#7A00B8", "#D47FFF", "#4C6BE3", "#A8BAFF", "#AAAAAA")

pdf(file='DataFigures/Fig6_FS2/cells_remap_prop.pdf', 2, 4, useD=F)
boxplot(percent_remap, ylim=c(0, 1), axes=F, xlim=c(0.6, 1.6), main='remapping cells', ylab='proportion of reliable cells')
points(rep(1,8), percent_remap, pch=21, bg=ratcols, cex=1)
axis(2, las=2)
dev.off()

#######################################################################
## normalised magnitude of change
cells_remap_reli <- cells_remap[(cells_remap$significance > -1)&(cells_remap$rate<150),]
p <-ggplot(cells_remap_reli, aes(x=rat.day, y=cross_norm.difference, col=rat.day)) + geom_violin(trim=TRUE)		
p + theme_minimal() + geom_jitter(shape=16, position=position_jitter(0.1)) + scale_color_manual(values=ratcols)

cells_remap_sign <- cells_remap[(cells_remap$significance > 0)&(cells_remap$rate<150),]

pdf(file='DataFigures/Fig6_FS2/cells_remap_sign.pdf', 6, 4, useD=F)
p <-ggplot(cells_remap_sign, aes(x=rat.day, y=cross_norm.difference, col=rat.day)) + geom_violin(trim=TRUE)		
p + theme_minimal() + geom_jitter(shape=16, position=position_jitter(0.1)) + scale_color_manual(values=ratcols)
dev.off()


