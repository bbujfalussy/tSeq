
dir.create('./DataFigures/Fig2/', showW=F)
chunkfile <- 'r1d2_20ms_decode.RData'
if (chunkfile %in% list.files()){
	load(file=chunkfile)
} else {
	###############################################
	### calculate the 20ms binned rasters with 5ms shift in each run episode
	chunks <- get_20ms_seqs(sp=spt, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_bins=theta_start, cellids=act_Ecells)
	
	###############################################
	### decode the bins
	Nchunks <- length(chunks)
	
	for (i_chunk in 1:Nchunks){
		bin_width <- chunks[[i_chunk]]$pos[,2] - chunks[[i_chunk]]$pos[,1]
		estimpos <- decode_ratemaps(chunks[[i_chunk]]$rast, bin_width, ratemaps, posmap)
		chunks[[i_chunk]]$estimpos <- estimpos
		cat(i_chunk, ' ')	
	}
	
	save(chunks, file=chunkfile)
}

###############################################
## for each point measure the relative error, theta phase and number of spikes
require(colormap)
require(ellipse)
require(abind)
require(viridis)


refx <- chunks[[1]]$pos[,4]
refy <- chunks[[1]]$pos[,5]
estx <- chunks[[1]]$estimpos$mean[,1]
esty <- chunks[[1]]$estimpos$mean[,2]
ref_dir <- chunks[[1]]$pos[,6]

errxy <- sqrt((refx-estx)^2 + (refy-esty)^2)
abs_estim_dir <- atan2(esty - refy, estx - refx)
estim_dir <- (abs_estim_dir - ref_dir + pi) %% (2*pi) - pi										
xi <- errxy * cos(estim_dir)
yi <- errxy * sin(estim_dir)

phase <- chunks[[1]]$pos[,7] %% (2*pi)
# plot(ncycles, col=viridis(100)[floor(phase*100)+1])
nsp <- rowSums(chunks[[1]]$rast, na.rm=T)

post_sd <- sqrt(rowSums(chunks[[1]]$estimpos$var))

aligned_xypns <- cbind(xi, yi, phase, nsp, post_sd)

Nchunks <- length(chunks)
for (i_chunk in 2:Nchunks){
	refx <- chunks[[i_chunk]]$pos[,4]
	refy <- chunks[[i_chunk]]$pos[,5]
	estx <- chunks[[i_chunk]]$estimpos$mean[,1]
	esty <- chunks[[i_chunk]]$estimpos$mean[,2]
	ref_dir <- chunks[[i_chunk]]$pos[,6]
	
	errxy <- sqrt((refx-estx)^2 + (refy-esty)^2)
	abs_estim_dir <- atan2(esty - refy, estx - refx)
	estim_dir <- (abs_estim_dir - ref_dir + pi) %% (2*pi) - pi										
	xi <- errxy * cos(estim_dir)
	yi <- errxy * sin(estim_dir)
	
	phase <- chunks[[i_chunk]]$pos[,7] %% (2*pi)
	# plot(ncycles, col=viridis(100)[floor(phase*100)+1])
	nsp <- rowSums(chunks[[i_chunk]]$rast, na.rm=T)
	
	post_sd <- sqrt(rowSums(chunks[[i_chunk]]$estimpos$var))

	aligned_xypns <- rbind(aligned_xypns, cbind(xi, yi, phase, nsp, post_sd))
}


## simple checks - similar counts at each phase, posterior variance decreases with spike count

M <- 20
brs <- 0:M/M*2*pi
hist(aligned_xypns[,3], brs, xlab='phase', main='')

hist(aligned_xypns[,5], 20, xlab='posterior SD', main='')
plot(aligned_xypns[,4], aligned_xypns[,5], pch='.', xlab='spike count', ylab='posterior SD')
plot(aligned_xypns[,3], aligned_xypns[,4], pch='.', xlab='phase', ylab='spike count')


###################################################
## the spike counts are different for each phase
## we perform thinning - downsampling data to match the spike counts distributions
thin_vec <- function(nsp, counts_thin){
	# nsp: spike counts in each bin to downsample
	# counts_thin: histogram, the number of elements to keep from 0 to Nmax
	Nmax <- length(counts_thin) - 1
	M <- length(counts_thin)
	index_keep <- NA
	for (j in 1:M){
		if (counts_thin[j] > 0){
			index_j <- which(nsp==(j-1))
			index_keep <- c(index_keep, sample(index_j, counts_thin[j]))
		}
	}
	sort(index_keep)
}

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

hnsp <- matrix(NA, M, 46)
for (i in 1:M){
	i_phase <- ((aligned_xypns[,3] >= brs[i]) & aligned_xypns[,3] < brs[i+1])
	aligned_xypns_i <- aligned_xypns[i_phase,]
	hnsp[i,] <- hist(aligned_xypns_i[,4], br=0:46-0.5, plot=F)$counts
}

image(0:19*18, 0:45, hnsp, axes=F, xlab='phase (deg)', ylab='spikes')
axis(2, las=2, 0:5*9)
axis(1, 0:4*90)
lines(0:19*18, hnsp %*% 0:45 / rowSums(hnsp), t='l')


#########################################################
## analyse decoded representation on the equalised histograms
th <- quantile(aligned_xypns[,4], 0.9)

rowSums(hnsp)
counts_thin <- apply(hnsp, 2, min)
thinning <- T

ms <- matrix(NA, 2, M)
sds <- matrix(NA, 2, M)
nsp <- matrix(NA, 2, M)
m_post_sd <- rep(NA, M)
xy_list <- list()


for (i in 1:M){
	i_phase <- ((aligned_xypns[,3] >= brs[i]) & aligned_xypns[,3] < brs[i+1])
	aligned_xypns_i <- aligned_xypns[i_phase,]
	# th <- quantile(aligned_xypns_i[,4], qq)
	if (thinning == T) {
		j_thin <- thin_vec(nsp=aligned_xypns_i[,4], counts_thin)
		aligned_xypns_i <- aligned_xypns_i[j_thin,]
	}	
	i_hsc <- aligned_xypns_i[,4] >= th
	xy_list[[i]] <- aligned_xypns_i[i_hsc,1:2]
	ms[,i] <- colMeans(aligned_xypns_i[i_hsc,1:2])
	sds[,i] <- apply(aligned_xypns_i[i_hsc,1:2], 2, sd)
	nsp[,i] <- c(th, mean(aligned_xypns_i[i_hsc,4]))
	hnsp[i,] <- hist(aligned_xypns_i[,4], br=0:46-0.5, plot=F)$counts
	m_post_sd[i] <- mean(aligned_xypns_i[i_hsc,5])
}

# image(0:19*18, 0:45, hnsp, axes=F, xlab='phase (deg)', ylab='spikes')
# axis(2, las=2, 0:5*9)
# axis(1, 0:4*90)
# lines(0:19*18, hnsp %*% 0:45 / rowSums(hnsp), t='l')
spread <- rep(NA, 20)
for (ii in 1:20) spread[ii] <- det(cov(xy_list[[ii]]))^(1/4)


pdf(paste('./DataFigures/Fig2/tseq20ms_bias_spread_th', th, '_rat', i_rat, '_day', i_day, '.pdf', sep=''), 6, 3, useD=F)
fillcols <- shifter(colormap(colormap=colormaps$hsv, 20, alpha=0.1), -7)
bgcols <- shifter(colormap(colormap=colormaps$hsv, 20), -7)
par(mfcol=c(1,2))
par(mar=c(4,4,4,4))
plot(0:19*18, ms[1,], xlim=c(0, 360), axes=F, pch=21, bg=bgcols, t='o', xlab='phase (deg)', ylab='error (cm)', ylim=c(-5, 20))
points(0:19*18, ms[2,], pch=22, bg=bgcols, t='o')
axis(1, c(0, 180, 360))
axis(2, las=2)

points(0:19*18, spread, pch=21, bg=grey(0.65), t='o')
# points(0:19*18, sqrt(sds[1,]^2 + sds[2,]^2)-20, pch=21, bg=2, t='o')
# points(0:19*18, sds[1,]-20, pch=21, bg=viridis(5)[4], t='o')
# points(0:19*18, sds[2,]-20, pch=21, bg=viridis(5)[5], t='o')
# points(0:19*18, nsp[2,], t='l', lwd=2)
# points(0:19*18, m_post_sd, t='l', lwd=2, col=3)
ell2 <- ellipse::ellipse

par(mar=c(2,2,2,2))
plot(xy_list[[6]][,1], xy_list[[6]][,2], axes=F, xlim=c(-20, 40), ylim=c(-30, 30), col=fillcols[6], pch=16, xlab='', ylab='')
for (i in c(8,10, 12,14, 16)){
	points(xy_list[[i]][,1], xy_list[[i]][,2], col=fillcols[i], pch=16)
}
for (i in c(6,8,10, 12,14, 16)){
	lines(ell2(cov(xy_list[[i]]), centre=colMeans(xy_list[[i]]), level=0.5), col=bgcols[i], t='l', lwd=2)
}
scalebar2(5, 5, '5 cm')
points(0, 0, pch=3)
dev.off()

######################################################
### real trajectory divergence
######################################################

NN <- dim(chunks[[2]]$pos)[1]
k <- seq(1, NN, by=20) # position by 100 ms

refx <- chunks[[2]]$pos[k,4]
refy <- chunks[[2]]$pos[k,5]
ref_dir <- chunks[[2]]$pos[k,6]

L <- length(k)
aligned_trajs <- array(NA, dim=c(L-20, 2, 21))
for (t_index in 11:(L-10)){
	ii <- (t_index-10):(t_index+10)
	refxy <- rbind(refx[ii] - refx[t_index], refy[ii] - refy[t_index])
	alpha <- (-1) * ref_dir[t_index]
	rotM <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), 2)
	aligned_trajs[t_index-10,,] <- rotM %*% refxy
	if (t_index==11){
		plot(aligned_trajs[t_index-10,1,], aligned_trajs[t_index-10,2,], col=viridis(21), t='o', xlab='xcoord', ylab='ycoord', xlim=c(-30, 30), ylim=c(-30, 30))		
	} else {
		lines(aligned_trajs[t_index-10,1,], aligned_trajs[t_index-10,2,], col=viridis(21), t='o')
	}
}

for (i_chunk in 3:Nchunks){
	NN <- dim(chunks[[i_chunk]]$pos)[1]
	if (NN > 600){
		k <- seq(1, NN, by=20) # position by 100 ms
		
		refx <- chunks[[i_chunk]]$pos[k,4]
		refy <- chunks[[i_chunk]]$pos[k,5]
		ref_dir <- chunks[[i_chunk]]$pos[k,6]
		
		L <- length(k)
		aligned_ts <- array(NA, dim=c(L-20, 2, 21))
		for (t_index in 11:(L-10)){
			ii <- (t_index-10):(t_index+10)
			refxy <- rbind(refx[ii] - refx[t_index], refy[ii] - refy[t_index])
			alpha <- (-1) * ref_dir[t_index]
			rotM <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), 2)
			aligned_ts[t_index-10,,] <- rotM %*% refxy
			# plot(refxy[1,], refxy[2,], col=viridis(21), xlim=c(-30, 30), ylim=c(-30, 30))
			# points(aligned_ts[t_index-10,1,], aligned_ts[t_index-10,2,], col=viridis(21))
		}
		aligned_trajs <- abind(aligned_trajs, aligned_ts, along=1)		
	}
}

ms_traj <- apply(aligned_trajs, c(2,3), mean)
sds_traj <- apply(aligned_trajs, c(2,3), sd)
spread_traj <- rep(NA, 21)
for (ii in 1:21) spread_traj[ii] <- det(cov(aligned_trajs[,,ii]))^(1/4)

pdf(paste('./DataFigures/Fig2/trajectory_divergence_rat', i_rat, '_day', i_day, '.pdf', sep=''), 6, 3, useD=F)
par(mfcol=c(1,2))
par(mar=c(4,4,4,4))
plot(-10:10/10, ms_traj[1,], xlim=c(-1, 1), axes=F, pch=21, bg=viridis(21), t='o', xlab='time (s)', ylab='distance (cm)', ylim=c(-20, 20))
points(-10:10/10, ms_traj[2,], pch=22, bg=viridis(21), t='o')
axis(1)
axis(2, las=2)

points(-10:10/10, spread_traj, pch=21, bg=grey(0.65), t='o')
legend('bottomright', c('forward', 'lateral', 'spread'), pch=c(21,22,21), pt.bg=c(viridis(4)[3], viridis(4)[3], grey(0.65)), bty='n', cex=0.8)

par(mar=c(2,2,2,2))
for (ii in seq(1, dim(aligned_trajs)[1], by=500)){
	if (ii ==1){
		plot(aligned_trajs[ii,1,], aligned_trajs[ii,2,], pch=16, col=viridis(21, alpha=0.2), t='o', xlab='', ylab='', xlim=c(-30, 30), ylim=c(-30, 30), axes=F)	
	} else {
		lines(aligned_trajs[ii,1,], aligned_trajs[ii,2,], pch=16, col=viridis(21, alpha=0.2), t='o')
	}
}
ii_plot <- 0:5*4+1
points(ms_traj[1, ii_plot], ms_traj[2, ii_plot], pch=21, bg=viridis(21)[ii_plot], cex=1)
for (ii in ii_plot){
	lines(ell2(cov(aligned_trajs[,,ii]), level=0.5, centre=ms_traj[,ii]), col=viridis(21)[ii], lwd=2)
}
scalebar2(5, 5, '5 cm')
points(0,0,pch=3, cex=2)
dev.off()


###############################################
## plot some decoded trajectories
if ((i_rat == 1) & (i_day == 2)){
	candidate_cycles <- 56:61
	i_chunk <- 144
	n_subplots <- length(candidate_cycles)
	tstarts <- rep(NA, n_subplots)
	pdf(file=paste('DataFigures/Fig2/R1D2_seq_', i_chunk, '_g2.pdf', sep=''), 3* n_subplots, 3, useD=F)
	par(mfcol=c(1,n_subplots))
	par(mar=c(1,1,1,1))
	chunk_rast <- list()
	
	for (i_sub in 1:n_subplots){
		i_cycle <- candidate_cycles[i_sub]
		estimpos <- chunks[[i_chunk]]$estimpos
		pos <- chunks[[i_chunk]]$pos
		kk <- which(pos[,3] == i_cycle)
		M <- length(kk)
		MM <- max(pos[,3])
		if (i_sub == 1) {
			tstart <- min(pos[kk,1])
			istart <- min(kk)
		}
		if (i_sub == n_subplots) {
			tend <- max(pos[kk,2])
			iend <- max(kk)			
		}
		tstarts[i_sub] <- min(pos[kk,1])
		chunk_rast[[i_sub]] <- chunks[[i_chunk]]$rast[kk,]
		# plot(pos[,4], pos[,5], xlim=c(0, 200), ylim=c(0, 200), pch=16, col=colormap(colormaps$grey, MM+2)[pos[,3]], axes=F, xlab='', ylab='')
		plot(pos[,4], pos[,5], xlim=c(0, 200), ylim=c(0, 200), t='l', col=grey(0.5), axes=F, xlab='', ylab='', lwd=1)
		box()
		estimSD <- sqrt(rowMeans(estimpos$var[kk,]))
		cex <- round((150 - estimSD) / 75, 1)
		lines(estimpos$mean[kk,1], estimpos$mean[kk,2], lwd=2, col=1)	
		points(estimpos$mean[kk,1], estimpos$mean[kk,2], pch=21, col=1, bg=grey(1), cex=cex)	
		# points(estimpos$mean[kk[3:(M-2)],1], estimpos$mean[kk[3:(M-2)],2], pch=21, col=1, bg=viridis(M-2, 1, 0.4, 1, option='B', dir=-1)[1:(M-2)], cex=cex[3:(M-2)])	
		points(estimpos$mean[kk[3:(M-2)],1], estimpos$mean[kk[3:(M-2)],2], pch=21, col=1, bg=gyr(M-4), cex=cex[3:(M-2)])	
		points(mean(pos[kk,4]), mean(pos[kk,5]), pch=3, col=viridis(10, option='D')[7], lwd=2, cex=2)
	}
	dev.off()
	
	jcells <- c(9, 58, 154, 176, 190, 261)
	
	png(file='DataFigures/Fig2/ratemaps6cells.png', 200, 1200)
	par(mfcol=c(6,1))
	par(mar=c(1,1,1,1)/2)
	for (j in jcells){
		image(attr(ratemaps, 'xcenters'), attr(ratemaps, 'ycenters'), ratemaps[j,,], axes=F, col=viridis(40, option='B'))
		points(pos[,4], pos[,5], pch=16, col=grey(0.75, alpha=0.25), t='o', cex=1)
	}
	dev.off()
	

	pdf(file='DataFigures/Fig2/R1D2_t33704_spikes.pdf', 9, 6, useD=F)
	layout(matrix(c(1,2), 2, 1), heights=c(2,1))
	par(mar=c(1,4.5, 1, 1))
	
	
	kk <- (spt[,1] > tstart) & (spt[,1] < tend)
	spikes_k <- spt[kk,]
	
	mm <- spikes_k[,2] %in% jcells
	spikes <- spikes_k[mm,]
	# plot(spikes_k[,1], spikes_k[,2], pch=16, ylim=c(0, 263), col=viridis(300, option='B', alpha=0.25)[spikes_k[,2]], xlab='', ylab='cells', axes=F, cex=0.5)
	# points(spikes[,1], spikes[,2], pch=16, ylim=c(0, 263), col=viridis(300, option='B', alpha=1)[spikes[,2]], lwd=1, cex=0.5)
	plot(spikes_k[,1], spikes_k[,2], pch='|', ylim=c(0, 263), col=viridis(300, option='D', alpha=1)[spikes_k[,2]], xlab='', ylab='cells', axes=F, cex=0.5)
	points(spikes[,1], spikes[,2], pch=15, ylim=c(0, 263), col=viridis(300, option='D', alpha=1)[spikes[,2]], lwd=1, cex=0.5)
	abline(v=tstarts)
	axis(2, las=2)
	
	par(mar=c(4.5,4.5, 1, 1))
	
	kk <- (lfp[,1] > tstart) & (lfp[,1] < tend)
	lfp_k <- lfp[kk,]
	plot(lfp_k[,1], lfp_k[,2], t='l', xlab='time (s)', ylab='LFP (uV)', axes=F)
	lines(lfp_k[,1], lfp_k[,4], t='l', col=2)
	axis(1)
	axis(2, las=2)
	
	dev.off()
}
