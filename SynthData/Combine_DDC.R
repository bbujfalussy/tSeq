########################################################################

library(viridis)
library(gplots)
require(colormap)


	
#######################################################################		
## decoding the spike rasters ...
codes <- c('MAP', 'DDC', 'sampling', 'PPC_var')

m_sigma <- array(NA, dim=c(4, 3), dimnames=list(codes, c('early', 'mid', 'late')))
se_sigma <- array(NA, dim=c(4, 3), dimnames=list(codes, c('early', 'mid', 'late')))

refcols <- c(rgb(95, 188, 211, max=255), rgb(109, 184, 1, max=255), rgb(254, 184, 1, max=255), rgb(153, 85, 255, max=255))
darkcols <- c(rgb(41, 127, 149, max=255), rgb(54, 92, 1, max=255), rgb(143, 103, 1, max=255), rgb(68, 0, 170, max=255))
lightcols <- c(rgb(193, 229, 238, max=255), rgb(181, 254, 75, max=255), rgb(255, 226, 153, max=255), rgb(218, 194, 255, max=255))

threshold <- '0' # 0 or median
Pvals <- matrix(NA, 4, 2, dimnames=list(codes, c('early-late', 'mid-late')))
Pfile <- 'DDC_Pvals.txt'
cat('P values for DDC decoding - one-sided, two-sample Kolmogorov-Smirnov test \n', file= Pfile)

for (threshold in c(0, 'median')){
	
	CDFfile <- paste('./SimFigs/Fig5/CDF_DDC_rate_th_', threshold, '_SDmin3cm_sim.pdf', sep='')
	cat('map used for decoding: rate, spike threshold: ', threshold, '\n', file= Pfile, append=T)
	
	pdf(file=CDFfile, 8, 3, useD=F)
	par(mfcol=c(1, 4))
	par(mar=c(4,5,3,2))
	for (i_code in 1:4){
			
		oparfile <- paste('./SimFigs/Fig5/decode_estR_Tmax1800_Ncells200_',codes[i_code], '_spPth32irreg_DDC.RData', sep='')
		load(oparfile)
		opars <- decode_DDC$opars
		###########################################################################
		## plot the data
			
		th.nsp <- 0 #median(opars[,4,]) # 0
		if (threshold == 'median') th.nsp <- median(opars[,4,])
	
		th_SD <- 200
		
		ii.gr1 <- opars[1,4,] > th.nsp
		ii.gr2 <- opars[2,4,] > th.nsp
		ii.gr3 <- opars[3,4,] > th.nsp
		
		CDF_SD_s <- cbind(cumsum(hist(opars[1,3,ii.gr1], breaks=3:201-0.5, plot=F)$density), cumsum(hist(opars[2,3,ii.gr2], breaks=3:201-0.5, plot=F)$density), cumsum(hist(opars[3,3,ii.gr3], breaks=3:201-0.5, plot=F)$density))
				
		ms_s_th <- rep(NA, 3)
		ss_s_th <- rep(NA, 3)
		
		j1 <- opars[1,3,] < th_SD
		j2 <- opars[2,3,] < th_SD
		j3 <- opars[3,3,] < th_SD
		Pvals[i_code, 1] <- ks.test(opars[1,3,j1& ii.gr1], opars[3,3,j3& ii.gr3], alternative='greater')$p
		Pvals[i_code, 2] <- ks.test(opars[2,3,j2& ii.gr2], opars[3,3,j3& ii.gr3], alternative='greater')$p
		cat(codes[i_code], 'early-late: ',  Pvals[i_code,1], 'mid-late: ',  Pvals[i_code,2], '\n', file= Pfile, append=T)
		
		
		{
			m_sigma[i_code, 1] <- mean(opars[1,3,j1& ii.gr1])
			se_sigma[i_code, 1] <- sd(opars[1,3,j1& ii.gr1]) / sqrt(sum(j1& ii.gr1))
			
			m_sigma[i_code, 2] <- mean(opars[2,3,j2& ii.gr2])
			se_sigma[i_code, 2] <- sd(opars[2,3,j2& ii.gr2]) / sqrt(sum(j2& ii.gr2))
			
			m_sigma[i_code, 3] <- mean(opars[3,3,j3& ii.gr3])
			se_sigma[i_code, 3] <- sd(opars[3,3,j3& ii.gr3]) / sqrt(sum(j3& ii.gr3))
		}
		
		
		matplot(3:200, CDF_SD_s, t='l', ylim=c(0.75, 1), lty=1, col=c(lightcols[i_code], refcols[i_code], darkcols[i_code]), lwd=c(4, 2,1 ), xlim=c(0, 50), xlab='', ylab='', main=paste(codes[i_code], ' - ', maps[calc_sharp+1], sep=''), axes=F, cex.lab=1.5, cex.main=1.5)
		mtext('cumulative probability', 2, cex=1, line=4)
		mtext('SD (cm)', 1, line=3, cex=1)
		axis(1, 0:2*25, cex.axis=1.5)
		axis(2, las=2, cex.axis=1.5)
		if (i_code == 1) {
			legend('bottomright', leg=c('early', 'mid', 'late'), lty=1, lwd=c(4,2,1), col=c(lightcols[1], refcols[1], darkcols[1]), bty='n', cex=1.5)
		}
	}
	cat( '\n', file= Pfile, append=T)
}


# dev.off()

cols <- matrix(c(lightcols, refcols, darkcols), 4)


means_file <- paste('./SimFigs/Fig5/means_DDC_rate_th_', threshold, '_SDmin3cm_sim.pdf', sep='')

ylims <- c(3, 12)
# if (calc_sharp == 1) ylims <- c(4, 8)
# if (calc_sharp == 2) ylims <- c(5, 10)

pdf(file=means_file, 2.5, 3, useD=F)
par(mar=c(4,4,1.6,1))
plot(1:3+1/10, m_sigma[1,], t='l', col=cols[1, 2], xlim=c(0.8, 3.8), ylim= ylims, axes=F, xlab='', ylab='decoded SD (cm)')
points(1:3+1/10, m_sigma[1,], pch=22, bg=cols[1,], cex=1.75)
plotCI(1:3+1/10, m_sigma[1,], se_sigma[1,], gap=0, pch=22, add=T, cex=1.75)

for (i_code in 2:4) {
	lines(1:3+ i_code/10, m_sigma[i_code,], t='l', col=cols[i_code, 2])
	points(1:3+ i_code/10, m_sigma[i_code,], pch=22, bg=cols[i_code,], cex=1.75)
	plotCI(1:3 + i_code/10, m_sigma[i_code,], se_sigma[i_code,], gap=0, pch=22, add=T, cex=1.75)

}

axis(1, 1:3+0.1, c('early', 'mid', 'late'), tick=F)
axis(2, 1:4*3, las=2)
dev.off()
