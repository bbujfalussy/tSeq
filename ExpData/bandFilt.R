bandFilt <- function(x, dt, f.low=NULL, f.high=NULL, graphics=F, fast=F){
	## x <- lfp.out; dt <- dt.out; graphics <- F
	## band-pass filter between frequencies f.low <= w <= f.high
	## the filter cutoff is a Gaussian cdf

	## 1. prepare a frequency axis for the FFT of the signal
	if (fast){
		LL <- length(x)
		L <- 2^ceiling(log(LL,2))
		xx <- x
		x <- rep(0, L)
		x[1:LL] <- xx
	} else {
		L <- length(x)	
	}
		
	Tmax <-  L * dt
	t <- seq(dt, Tmax, by=dt)
	dfreq <- 1/Tmax
	freq <- seq(0, (L-1)/L/dt, by=1/Tmax)
	n0 <- floor(L/2)
	n1 <- ceiling(L/2) + 2
	freq[n1:L] <- rev(freq[2:n0])

	## 2. fft of the signal
	fft.x <- fft(x)

	# NN <- 12/dfreq
	# plot(freq[1:NN],  log(abs(fft.x[1:NN])^2), t="h", xlab="frequency", ylab="log power")
	# points(freq[L:(L-NN+1)], log(abs(fft.x[L:(L-NN+1)])^2), col=2)

	## 3. the filter kernel
	if (is.null(f.low) & is.null(f.high)) stop("either f.low of f.high must be provided!")
	if (is.null(f.low)) f.low <- 0
	if (is.null(f.high)) f.high <- max(freq)
	if (f.low > f.high) stop("f.low must be smaller than f.high!")
	
	w.freq <- 3
	if (f.low > 4) w.freq <- max(w.freq, round(0.2*min(f.low, f.high-f.low)))
	if ((f.high - f.low) < 3) w.freq <- 2
	if ((f.high - f.low) < 2) w.freq <- 0

	ifilt <- rep(0, L)
	i.true <- (freq>=f.low)&(freq<=f.high)
	ifilt[i.true] <- 1

	ii.up.1 <- which((freq - (f.low - w.freq/2))^2 == min((freq - (f.low - w.freq/2))^2))[1]
	ii.up.2 <- which((freq - (f.low + w.freq/2))^2 == min((freq - (f.low + w.freq/2))^2))[1]
	ii1 <- ii.up.1:ii.up.2; LL1 <- length(ii1); nn1 <- seq(-3, 3, length=LL1)
	ifilt[ii1] <- pnorm(nn1)
	
	ii.down.1 <- which((freq - (f.high - w.freq/2))^2 == min((freq - (f.high - w.freq/2))^2))[1]
	ii.down.2 <- which((freq - (f.high + w.freq/2))^2 == min((freq - (f.high + w.freq/2))^2))[1]
	ii2 <- ii.down.1:ii.down.2; LL2 <- length(ii2); nn2 <- seq(-3, 3, length=LL2)
	ifilt[ii2] <- rev(pnorm(nn2))

	if (f.low == 0){
		ifilt[ii1] <- 1
		ifilt[L - ii1 + 1] <- 1
		ifilt[L-ii2+1] <- rev(pnorm(nn2))
	} else {
		ifilt[L-ii1+2] <- pnorm(nn1)
		ifilt[L-ii2+2] <- rev(pnorm(nn2))
	}

	# plot(freq, ifilt, xlab="frequency", ylab="filter power")


	# fft2.filt <- Re(fft(ifilt, inverse=T))/L
	# plot(t[1:200], fft2.filt[1:200], t="l")
	# plot(t, fft2.filt, t="l")

	## 4. the inverse FFT
	
	fft2.x <- Re(fft(fft.x*ifilt, inverse=T))/L
	if (fast){
		fft2.x <- fft2.x[1:LL]
	}

	if (graphics==T){
		plot(t[1:LL], x[1:LL], t="l")
		lines(t[1:LL], fft2.x, col=2)
	}
	
	fft2.x
}