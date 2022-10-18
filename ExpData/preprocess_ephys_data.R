########################################################################
### load position, spikes and LFP, filter the LFP, detect artefats and calculate casc ratemaps

cat('loading data for rat', i_rat, 'day', i_day, '... \n')

ddata <- data[[i_rat]][[i_day]]
ncells <- max(ddata$spikes[,2])

## filter position to remove noise and make it smoother
txyvd <- posfilter(ddata$position, sdfilt=0.25, graphics=T)

spt <- ddata$spikes

tmin <- min(txyvd[,1])
tmax <- max(txyvd[,1])

## detect periods of movement - at least 1s long move with at least 0.2 s stops between
dt_pos <- median(diff(txyvd[,1]))
vmin <- 5 # cm/s
T_run_min <- 1 # s
T_stop_min <- 0.2 # s
ii_run <- txyvd[,4] > vmin
ii_run[is.na(ii_run)] <- F

xyvd <- txyvd[,2:5]
xyvd[!ii_run] <- NA

dir.create('./DataFigures/Fig3_FS2/', showW=F)
motionfile <- paste('DataFigures/Fig3_FS2/R', i_rat, 'D', i_day, '_motion_stats.pdf', sep='') else file=NULL
stats <- plot_motion_stats(xyvd=xyvd, dt_pos= dt_pos, file=motionfile)

#####################################
## we choose the lfp where the theta is stronger - larger correlation between the raw and the theta filtered
## Rat					1		1		2		2		3		3		4		4
## Day					1		2		1		2		1		2		1		2
## Left of Right		L		L		R		L		L		L?		R		L-R
## number of cells	213	263	181	164	80		80		97		186
## theta				+		+		-		-		+		-		+		++
## artefacts			-		-		+		++		-		-		-		±   
cat('loading LFP data \n')
if ((10 * i_rat + i_day) %in%  c(11, 12, 22, 31, 32) ){
# if ((10 * i_rat + i_day) %in%  c(23, 41, 42) ){
	lfp_side <- 'Left'
	lfp_name <- 'LFP.Left'
} else {
	lfp_side <- 'Right'
	lfp_name <- 'LFP.Right'
}
lfp_file <- paste('./Pfeiffer13/Rat', i_rat, '_Day', i_day, '_LFP_', lfp_side, '.mat', sep='')
tryCatch({lfp_data <- readMat(lfp_file)}, error = function(e) print('error: LFP file not found. Please contact the authors to obtain the data.'), warning=function(w) print('error: LFP file not found.  Please contact the authors to obtain the data.'))
freq <- as.vector(lfp_data$LFP.Frequency)
dt <- 1/freq # s


################################################################
## extracting band-pass filtered lfp signals to detect theta and artefacts...
imin <- min(which(lfp_data[[lfp_name]][,1] > tmin))
imax <- max(which(lfp_data[[lfp_name]][,1] < tmax))
t_theta <- lfp_data[[lfp_name]][imin:imax,1]
rawlfp <- lfp_data[[lfp_name]][imin:imax,2]

L_lfp <- length(rawlfp)
r_lfp <- rep(0, 2^ceiling(log(L_lfp, 2)))
r_lfp[1:L_lfp] <- rawlfp
lfp_low <- ffilter(r_lfp, 1/dt, to=12, rescale=T)
lfp_low <- lfp_low[1:L_lfp]
### theta and ripple band
lfp_theta <- bandFilt(rawlfp, dt, f.low=4, f.high=12, graphics=F, fast=T)
lfp_ripple <- bandFilt(rawlfp, dt, f.low=120, f.high=250, graphics=F, fast=T)

## high frequency artefacts
a_ripple <- env(lfp_ripple, freq, plot=F)
ripple_low <- filter(c(rep(0, 320), a_ripple), rep(1/321, 321))
ripple_low  <- ripple_low[!is.na(ripple_low)]

lfp <- cbind(t_theta, rawlfp, lfp_ripple, lfp_theta)


## theta phase - from Hilbert trafo
# a_theta <- env(lfp_theta, freq, plot=F)
h_theta <- ifreq(lfp_theta, freq, phase=T, plot=F)$p[,2]


#######################
# # ARTEFACTS - removing recording artefacts - we mark them as 'no run' periods and do not analyse those sections
#######################
## no spiking for more than 200 ms (in certain recordings the spikes were not detected for some periods of time)

ISI <- diff(spt[,1])
isp_first <- which(ISI>0.2)
if (length(isp_first > 0)){
	for (ii in isp_first){
		t1 <- spt[ii,1]
		t2 <- spt[ii+1,1]
		if ((t2 > txyvd[1,1]) & (t1 < tail(txyvd[,1],1))){
			i1 <- which.min((txyvd[,1] - t1)^2) - 1 
			i2 <- which.min((txyvd[,1] - t2)^2) + 1
			cat(i2-i1, ' ')
			ii_run[i1:i2] <- F
		}
	}
}

## t_theta is recorded at lower resolution
print('LFP resolution artefacts...')
ITI <- which(diff(t_theta) > 0.0004)
NN <- length(ITI)

if (NN > 1){ # there are at least two breaks
	ii <- 1	# start at the first one
	while (ii <= NN) {
		t1 <- t_theta[ITI[ii]]
		jj <- 1
		while((ii + jj < NN) & (ITI[ii + jj] == ITI[ii]+jj)){
			jj <- jj + 1
		}
		jj <- jj - 1
		t2 <- t_theta[ITI[ii+jj]+1]				
		j1 <- which.min((txyvd[,1] - t1)^2) - 1 
		j2 <- which.min((txyvd[,1] - t2)^2) + 1 
		cat(ii, ii+jj, j2-j1, '\n ')
		ii_run[j1:j2] <- F
		ii <- ii + jj + 1
	}	
} else if (NN == 1){
	t1 <- t_theta[ITI]
	t2 <- t_theta[ITI+1]			
	j1 <- which.min((txyvd[,1] - t1)^2) - 1 
	j2 <- which.min((txyvd[,1] - t2)^2) + 1 
	cat(j2-j1, ' ')
	ii_run[j1:j2] <- F
}

## low pass LFP is very high
lfp_amp <- abs(approx(t_theta , lfp_low, txyvd[,1])$y)
for (i in which(lfp_amp > 1000)) {
	j <- max(1, i-5)
	jj <- min(length(ii_run), i+5)
	ii_run[j:jj] <- F
}

## high-pass LFP is very high 
ripple_amp <- abs(approx(t_theta , ripple_low, txyvd[,1])$y)
for (i in which(ripple_amp > 70)) {
	j <- max(1, i-5)
	jj <- min(length(ii_run), i+5)
	ii_run[j:jj] <- F
}


## detect short breaks in theta - where the theta frequency is shorter than 12Hz or longer than 5Hz
i_run <- get_long_periods(ii_run, T_run_min / dt_pos, T_stop_min / dt_pos)
print(sum(i_run) / length(i_run))


# tstart <- 2770
# pdf (file='R1D1_t2774_Left.pdf', 10, 5, useD=F)
# iplot.run(txyvd, lfp, spt, tmin=tstart, tmax=tstart+5)
# dev.off()		

##################################################################################
## estimate the ratemaps
########################################################################
cat('calculating ratemaps \n')

xy_brs <- seq(0, 200, by=5)

posmap <- estimate_posmap(txyvd[,1:3], i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, graphics=T)

rmaps <- array(NA, dim=c(ncells, 40, 40))
nspikes <- rep(NA, ncells)

par(mfcol=c(3,4))
for (ii in 1:ncells){
	spi <- spt[spt[,2] == ii,1]
	ratemap <- estimate_ratemap(spi, txyvd[,1:3], i_data=i_run, xbrs= xy_brs, ybrs=xy_brs, sigma=10, posmap=posmap, graphics=0)
	rmaps[ii,,] <- ratemap
	nspikes[ii] <- attr(ratemap, 'nsp')
}

rmaps[rmaps < 0.1] <- 0.1
act_cells <- which(nspikes > 200) 
# act_Ecells <- which(which(act_cells) %in% ddata$E_Cells)
act_Ecells <- ddata$E_Cells[ddata$E_Cells[,1] %in% act_cells,1]
ratemaps <- rmaps#[act_cells,,]

attr(ratemaps, 'xcenters') <- attr(posmap, 'xcenter')
attr(ratemaps, 'ycenters') <- attr(posmap, 'xcenter')

dir.create('./DataFigures/Fig3_FS3/', showW=F)

filename <- paste('DataFigures/Fig3_FS3/R', i_rat, 'D', i_day, '_place_stats.png', sep='')
png(file=filename, 1200, 600, pointsize=16)
FLB <- plot_rate_stats(ratemaps, act_Ecells, act_Ecells, c(9, 10, 18, 19), posmap=posmap) # cm, Fisher lower bound on decoding error
dev.off()
