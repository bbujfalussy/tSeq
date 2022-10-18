# chunks <- run_chunks(sp=ddata$spikes, txy=txyv[,1:3], i_data=i_run, dt=0.1, cellids=act_Ecells)
qq <- 19/20 # 19/20
theta_phases <- seq(0, by=pi/4, length=8)
decfilename <- paste('R', i_rat, 'D', i_day, '_q', qq, '_dec_chunks.RData', sep='')
if (!('Fig4' %in% list.files('DataFigures'))){
	dir.create('./DataFigures/Fig4/', showW=F)
}
if (decfilename %in% list.files('DataFigures/Fig4/')){
	cat('loading decoded theta chunks...')
	decfilename <- paste('./DataFigures/Fig4/', decfilename, sep='')
	load(decfilename)
} else {
	decfilename <- paste('./DataFigures/Fig4/', decfilename, sep='')
	cat('calculating theta chunks...')
	
	chunks1 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi =0, L_phi=2*pi/3, cellids=act_Ecells)
	chunks2 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi=pi/4, L_phi=2*pi/3, cellids=act_Ecells)
	chunks3 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi=2*pi/4, L_phi=2*pi/3, cellids=act_Ecells)
	chunks4 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi=3*pi/4, L_phi=2*pi/3, cellids=act_Ecells)
	chunks5 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi=pi, L_phi=2*pi/3, cellids=act_Ecells)
	chunks6 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi=5*pi/4, L_phi=2*pi/3, cellids=act_Ecells)
	chunks7 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi=6*pi/4, L_phi=2*pi/3, cellids=act_Ecells)
	chunks8 <- get_phase_chunks(sp=ddata$spikes, txy=txyvd,  i_data=i_run, phi=h_theta, t_theta=t_theta, start_phi=7*pi/4, L_phi=2*pi/3, cellids=act_Ecells)	
	
	cat('decoding position ...')
	estimpos1 <- decode_theta_chunks(chunks1, ratemaps, posmap, calc_error=T)
	estimpos2 <- decode_theta_chunks(chunks2, ratemaps, posmap, calc_error=T)
	estimpos3 <- decode_theta_chunks(chunks3, ratemaps, posmap, calc_error=T)
	estimpos4 <- decode_theta_chunks(chunks4, ratemaps, posmap, calc_error=T)
	estimpos5 <- decode_theta_chunks(chunks5, ratemaps, posmap, calc_error=T)
	estimpos6 <- decode_theta_chunks(chunks6, ratemaps, posmap, calc_error=T)
	estimpos7 <- decode_theta_chunks(chunks7, ratemaps, posmap, calc_error=T)
	estimpos8 <- decode_theta_chunks(chunks8, ratemaps, posmap, calc_error=T)
	
	dec_chunks <- list(chunks1=chunks1, chunks2=chunks2, chunks3=chunks3, chunks4=chunks4, chunks5=chunks5, chunks6=chunks6, chunks7=chunks7, chunks8=chunks8)
	dec_chunks$est1 <- estimpos1
	dec_chunks$est2 <- estimpos2
	dec_chunks$est3 <- estimpos3
	dec_chunks$est4 <- estimpos4
	dec_chunks$est5 <- estimpos5
	dec_chunks$est6 <- estimpos6
	dec_chunks$est7 <- estimpos7
	dec_chunks$est8 <- estimpos8
	save(dec_chunks, file=decfilename)
}
##################################################
## aligned decoded positions, number of spikes and bin duration

aligned_errs <- list()
rates_mean_sd <- matrix(NA, 8, 2)
n_cells <- sum(!is.na(dec_chunks$chunks1[[1]]$rast[1,]))
max_nsp <- 0

for (j in 1:8){
	xynt <- c(0, 0, 0, 0)
	if (j == 1) {estimpos <- dec_chunks$est1; inchunk <- dec_chunks$chunks1}
	if (j == 2) {estimpos <- dec_chunks$est2; inchunk <- dec_chunks$chunks2}
	if (j == 3) {estimpos <- dec_chunks$est3; inchunk <- dec_chunks$chunks3}
	if (j == 4) {estimpos <- dec_chunks$est4; inchunk <- dec_chunks$chunks4}
	if (j == 5) {estimpos <- dec_chunks$est5; inchunk <- dec_chunks$chunks5}
	if (j == 6) {estimpos <- dec_chunks$est6; inchunk <- dec_chunks$chunks6}
	if (j == 7) {estimpos <- dec_chunks$est7; inchunk <- dec_chunks$chunks7}
	if (j == 8) {estimpos <- dec_chunks$est8; inchunk <- dec_chunks$chunks8}
	

	for (i in 1:length(estimpos)){
		if (length(estimpos[[i]]$dir) > 1){
			chunk <- estimpos[[i]]
			xi <- chunk$rsqerror * cos(chunk$dir)
			yi <- chunk$rsqerror * sin(chunk$dir)
			# points(xi, yi, pch=16, cex=1, col=cols[j])
			nsp <- rowSums(inchunk[[i]]$rast, na.rm=T)
			t_bin <- inchunk[[i]]$pos[,5]
			xynt <- rbind(xynt, cbind(xi, yi, nsp, t_bin))
		}
	}
	if (max(xynt[,3]) > max_nsp) max_nsp <- max(xynt[,3])
	aligned_errs[[j]] <- xynt[-1,]
}


##################################################
## plot the decoded positions and thear mean and spread
qq <- 0/20

brs <- 0 :(max_nsp+1) - 0.5
hnsp <- matrix(NA, 8, max_nsp+1)
for (j in 1:8){
	hnsp[j,] <- hist(aligned_errs[[j]][,3], br=brs, plot=F)$counts
}
image(1:8, 0:max_nsp, hnsp, axes=F, ylab='', xlab='')
lines(1:8, hnsp %*% 0: max_nsp / rowSums(hnsp), lwd=2)
axis(2, las=2)
axis(1, 1:8, 1:8*45)
min_counts <- apply(hnsp, 2, min)
plot(0:max_nsp, cumsum(min_counts) / sum(min_counts))
spike_threshold <- min(which(cumsum(min_counts) / sum(min_counts) > qq)) - 1
abline(h=qq, v=spike_threshold)

## function to perform thinning
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

#############################################################################################
## to make a fair comparison, we need to perform thinning - downsample nspike-histogram to equaliye counts accross phase

cols <- colormap(colormap=colormaps$hsv, 8, alpha=0.25)
mcols <- colormap(colormap=colormaps$hsv, 18)[c(1, 3, 4, 5, 9, 11, 13, 15)]

ell2 <- ellipse::ellipse

ms <- matrix(NA, 8, 2)
spread <- rep(NA, 8)
rate <- matrix(NA, 8, 2)

spike_threshold <- min(which(cumsum(min_counts) / sum(min_counts) > qq)) - 1

filename <- paste('DataFigures/Fig4/R', i_rat, 'D', i_day, '_q', qq, '_postmean_phase.png', sep='')
if (graphics) png(file=filename, 1200, 1200, pointsize=48)

# filename <- paste('DataFigures/Fig4/R', i_rat, 'D', i_day, '_q', qq, '_postmean_phase.pdf', sep='')
# if (graphics) pdf(file=filename, 4, 4, useD=F)
par(mar=c(1,1,1,1))
xylist <- list()

plot(0, 0, axes=F, xlim=c(-20, 40), ylim=c(-30, 30), xlab='', ylab='', pch=3, cex=3)
for (j in 1:8){
	xynt <- aligned_errs[[j]]
	all_nsp <- xynt[,3]
	all_t_bin <- xynt[,4]
	jj <- thin_vec(all_nsp, min_counts)
	xy <- xynt[jj,1:2]
	nsp <- xynt[jj,3]
	ii <- which(nsp > spike_threshold)
	xylist[[j]] <- xy[ii,]
	ms[j,] <- colMeans(xy[ii,])
	spread[j] <- det(cov(xy[ii,]))^(1/4)
	rate[j,1] <- sum(all_nsp) / sum(all_t_bin) / n_cells	
	rate[j,2] <- sd(all_nsp / all_t_bin / n_cells)
	points(xy[ii,1], xy[ii,2], pch=16, cex=1, col=cols[j])
}
for (j in 1:8){
	lines(ell2(cov(xylist[[j]]), centre=c(ms[j,1], ms[j,2]), level=0.5), t='l', col=mcols[j], lwd=2)
	points(ms[j,1], ms[j,2], pch=21, col=1, bg=mcols[j], cex=2)	
}
points(0, 0, pch=3, cex=3)
scalebar2(5, 5, ' ', '5 cm')
if (graphics) dev.off()
	
			
initpars <- c(6, 6, 0)
errcos(initpars, theta_phases, ms[,1])
newpars <- optim(initpars, errcos, gr=NULL, xx=theta_phases, yy=ms[,1])
xx <- seq(0, 4*pi, length=100)
yy <- predcos(newpars$par, xx)
yy_ref <- predcos(c(1, 1, pi), xx)

filename <- paste('DataFigures/Fig4/R', i_rat, 'D', i_day, '_q', qq, '_theta_bias.pdf', sep='')
if (graphics) pdf(file=filename, useD=F, 3, 4)
par(mfcol=c(3,1))
par(mar=c(1, 4, 4, 1))

captiontext <- paste('peak of phase:', round(newpars$par[3], 2), 'rad')
plot(c(theta_phases, 2*pi+ theta_phases), c(ms[,1], ms[,1]), col=1, pch=21, t='o', bg=mcols, xlab='', ylab='decoding bias (cm)', axes=F, ylim=(ceiling(range(ms)) + c(-1, 0)), xlim=c(0, 4*pi), main=captiontext)
points(c(theta_phases, 2*pi + theta_phases), c(ms[,2], ms[,2]), col=1, pch=22, bg=mcols, t='o')

lines(xx, yy)
axis(2, las=2)
legend('topright', leg=c('forward', 'lateral'), pch=c(21, 22), pt.bg=mcols[2], bty='n')

plot(c(theta_phases, 2*pi + theta_phases), c(spread, spread), col=1, pch=21, t='o', bg=mcols, xlab='', ylab='spread (cm)', axes=F, xlim=c(0, 4*pi), ylim=range(spread), main='120 deg window')
axis(2, las=2)


par(mar=c(4, 4, 1, 1))
plotCI(c(theta_phases, 2*pi + theta_phases), c(rate[,1], rate[,1]), uiw=c(rate[,2], rate[,2]), gap=0, ylim=c(0, 3), axes=F, xlab='theta phase (deg)', ylab='firing rate (Hz)', xlim=c(0, 4*pi))
axis(1, c(0, pi, 2*pi, 3*pi, 4*pi), c(0, 180, 360, 180, 360))
axis(2, 0:3, las=2)
lines(xx, yy/max(yy) * max(rate[,1]))
points(c(theta_phases, 2*pi + theta_phases), c(rate[,1], rate[,1]), pch=21, t='o', bg=mcols)

if (graphics) {
	dev.off()

	bias_rate <- list(ms=ms, spread=spread, rate=rate)
	biasfilename <- paste('DataFigures/Fig4/R', i_rat, 'D', i_day, '_q', qq, '_theta_bias.RData', sep='')
	save(bias_rate, file=biasfilename)
}
