

###################################################################
### MAP decoding the spikes in the three different window of the theta
###################################################################

outfile <- paste('R', i_rat, 'D', i_day, '_decoded_theta_seqs.RData', sep='')
if (outfile %in% list.files()) {
	print('loading decoded position')
	load(file=outfile)
	chunks <- dec_theta$chunks
	estimpos_early <- dec_theta$est_early
	estimpos_mid <- dec_theta$est_mid
	estimpos_late <- dec_theta$est_late
} else {
## early mid and late theta phases... - for testing PPC, segments of similar spike numbers
	print('decoding position')
	chunks <- get_theta_seqs(sp=spt, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_bins=theta_start, cellids=act_Ecells)
	cat('chunks calculated... ')
	estimpos_early <- decode_theta_chunks(chunks$past, ratemaps, posmap, calc_error=T, xshift = 0, tshift = 0, txyv=txyvd)
	cat('decoded: past, ')
	estimpos_mid <- decode_theta_chunks(chunks$present, ratemaps, posmap, calc_error=T, xshift = 0, tshift = 0, txyv=txyvd)
	cat('present, ')
	estimpos_late <- decode_theta_chunks(chunks$future, ratemaps, posmap, calc_error=T, xshift = 0, tshift = 0, txyv=txyvd)
	cat('future. \n')
	dec_theta <- list(chunks=chunks, est_early=estimpos_early, est_mid=estimpos_mid, est_late=estimpos_late)
	save(dec_theta, file=outfile)
}

#########################################################
## rmaps based on early phase spikes with temporally shifted position data
#########################################################
mapsname <- paste('3maps_rat', i_rat, '_day', i_day, '.RData', sep='')
if (mapsname %in% list.files()) {
	print('loading ratemaps')
	load(file=mapsname)
	ratemaps <- maps$ratemaps
	ratemaps_eml <- maps$ratemaps_eml
} else{
	print('calculating ratemaps ... ')
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
	
	n3spikes <- colSums(rast_early, na.rm=T) + colSums(rast_mid, na.rm=T) + colSums(rast_late, na.rm=T)
	txy_early <- txy_early[!is.na(txy_early[,1]),]
	txy_mid <- txy_mid[!is.na(txy_mid[,1]),]
	txy_late <- txy_late[!is.na(txy_late[,1]),]

	rates_eml <- cbind(colSums(rast_early, na.rm=T) / sum(txy_early[,5]), colSums(rast_mid, na.rm=T) / sum(txy_mid[,5]), colSums(rast_late, na.rm=T) / sum(txy_late[,5]), nspikes / (sum(txy_late[,5]) + sum(txy_mid[,5]) + sum(txy_early[,5])))

	# use temporally shifted the position data to minimise error for early mid and late phase spikes
	dir.create('./DataFigures/Fig4_FS2/', showW=F)
	shiftfile <- paste('DataFigures/Fig4_FS2/R', i_rat, 'D', i_day, '_q', quantile_curvature, 'decoding_curves_eml.pdf', sep='')
	source('./decoding_curves_exp.R')
	
	xx_shift <- approx(txy_early[,1], txy_early[,2], xout= txy_early[,1] + delay_eml[1], rule=2)$y 
	yy_shift <- approx(txy_early[,1], txy_early[,3], xout= txy_early[,1] + delay_eml[1], rule=2)$y 
	xyw_early <- cbind(xx_shift, yy_shift, txy_early[,5])

	xx_shift <- approx(txy_mid[,1], txy_mid[,2], xout= txy_mid[,1] + delay_eml[2], rule=2)$y 
	yy_shift <- approx(txy_mid[,1], txy_mid[,3], xout= txy_mid[,1] + delay_eml[2], rule=2)$y 
	xyw_mid <- cbind(xx_shift, yy_shift, txy_mid[,5])

	xx_shift <- approx(txy_late[,1], txy_late[,2], xout= txy_late[,1] + delay_eml[3], rule=2)$y 
	yy_shift <- approx(txy_late[,1], txy_late[,3], xout= txy_late[,1] + delay_eml[3], rule=2)$y 
	xyw_late <- cbind(xx_shift, yy_shift, txy_late[,5])

	posmap_early <- estimate_posmap_xyw(xyw_early, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F, prior=1/100)
	posmap_mid <- estimate_posmap_xyw(xyw_mid, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F, prior=1/100)
	posmap_late <- estimate_posmap_xyw(xyw_late, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=F, prior=1/100)

	ratemaps_eml <- array(NA, dim=c(ncells, 3, 40, 40))
	Skaggs_eml <- array(NA, dim=c(ncells, 4))
	PFsize_eml <- array(NA, dim=c(ncells, 4))

	print(' ... from early, mid and late spikes ...')
	
	for (ii in 1:ncells){
		counts <- rast_early[,ii]
		ratemap <- estimate_ratemap_xyw(counts, xyw_early, xbrs=xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap_early, graphics=F)
		ratemaps_eml[ii,1,,] <- ratemap
		Skaggs_eml[ii, 1] <- skaggs93.info(ratemap, posmap_early) / rates_eml[ii, 1]
		PFsize_eml[ii, 1] <- sum(ratemap > 0.1*max(ratemap)) / 400
		
		counts <- rast_mid[,ii]
		ratemap <- estimate_ratemap_xyw(counts, xyw_mid, xbrs=xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap_mid, graphics=F)
		ratemaps_eml[ii,2,,] <- ratemap
		Skaggs_eml[ii, 2] <- skaggs93.info(ratemap, posmap_mid) / rates_eml[ii, 2]
		PFsize_eml[ii, 2] <- sum(ratemap > 0.1*max(ratemap)) / 400
	
		counts <- rast_late[,ii]
		ratemap <- estimate_ratemap_xyw(counts, xyw_late, xbrs=xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap_late, graphics=F)
		ratemaps_eml[ii,3,,] <- ratemap
		Skaggs_eml[ii, 3] <- skaggs93.info(ratemap, posmap_late) / rates_eml[ii, 3]
		PFsize_eml[ii, 3] <- sum(ratemap > 0.1*max(ratemap)) / 400
	
		Skaggs_eml[ii, 4] <- skaggs93.info(ratemaps[ii,,], posmap) / rates_eml[ii, 4]
		PFsize_eml[ii, 4] <- sum(ratemaps[ii,,] > 0.1*max(ratemaps[ii,,])) / 400
	}

	#########################################################
	## plot the map features baed on early, mid and late spikes
	#########################################################	
	
	dir.create('./DataFigures/Fig2/', showW=F)
	PFsize_file <- paste('DataFigures/Fig2/PFsize_R', i_rat, 'D', i_day, '.pdf', sep='')
	pdf(file=PFsize_file, 7, 4, useD=F)

	par(mfcol=c(1,2))
	plot(Skaggs_eml[act_Ecells, 1], Skaggs_eml[act_Ecells, 4], pch=16, cex=1, col=viridis(5, alpha=0.15)[3], xlab='early', ylab='mid or late', main='Skaggs info (bit/spike)')
	points(Skaggs_eml[act_Ecells, 1], Skaggs_eml[act_Ecells, 2], pch=16, cex=1, col=viridis(5, alpha=0.15,option='B')[3])
	points(Skaggs_eml[act_Ecells, 1], Skaggs_eml[act_Ecells, 3], pch=16, cex=1, col=viridis(5, alpha=0.15,option='B')[4])
	points(median(Skaggs_eml[act_Ecells, 1]), median(Skaggs_eml[act_Ecells, 4]), pch=3, cex=1.5, col=viridis(5,option='D')[3], lwd=2)
	points(median(Skaggs_eml[act_Ecells, 1]), median(Skaggs_eml[act_Ecells, 2]), pch=3, cex=1.5, col=viridis(5,option='B')[3], lwd=2)
	points(median(Skaggs_eml[act_Ecells, 1]), median(Skaggs_eml[act_Ecells, 3]), pch=3, cex=1.5, col=viridis(5,option='B')[4], lwd=2)
	abline(0,1, col=grey(0.75))
	
	plot(PFsize_eml[act_Ecells, 1], PFsize_eml[act_Ecells, 4], pch=16, cex=1, col=viridis(5, alpha=0.15)[3], xlab='early', ylab='mid or late', main='PF_size (m2)')
	points(PFsize_eml[act_Ecells, 1], PFsize_eml[act_Ecells, 2], pch=16, cex=1, col=viridis(5, alpha=0.15,option='B')[3])
	points(PFsize_eml[act_Ecells, 1], PFsize_eml[act_Ecells, 3], pch=16, cex=1, col=viridis(5, alpha=0.15,option='B')[4])
	points(median(PFsize_eml[act_Ecells, 1]), median(PFsize_eml[act_Ecells, 4]), pch=3, cex=1.5, col=viridis(5,option='D')[3], lwd=2)
	points(median(PFsize_eml[act_Ecells, 1]), median(PFsize_eml[act_Ecells, 2]), pch=3, cex=1.5, col=viridis(5,option='B')[3], lwd=2)
	points(median(PFsize_eml[act_Ecells, 1]), median(PFsize_eml[act_Ecells, 3]), pch=3, cex=1.5, col=viridis(5,option='B')[4], lwd=2)
	abline(0,1, col=grey(0.75))
	legend('bottomright', leg=c('mid', 'late', 'all'), pch=16, col=c(viridis(5, alpha=0.25,option='B')[3], viridis(5, alpha=0.25,option='B')[4], viridis(5, alpha=0.25,option='D')[3]), bty='n')
	
	dev.off()


	
	########################################################
	# saving the ratemaps ...
	########################################################
	maps <- list(ratemaps=ratemaps, ratemaps_eml= ratemaps_eml)
	maps$posmap <- posmap
	maps$posmap_early <- posmap_early
	maps$posmap_mid <- posmap_mid
	maps$posmap_late <- posmap_late
	maps$Skaggs_eml <- Skaggs_eml
	maps$PFsize_eml <- PFsize_eml
	maps$rates_eml <-rates_eml
	maps$act_Ecells <- act_Ecells
	save(maps, file=mapsname)
}


