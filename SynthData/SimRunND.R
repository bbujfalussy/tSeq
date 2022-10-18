### functions to generate synthetic population activity in the hippocampus
### firing of place cells while the rat explores a N-Dimensional environment
### 1. sim.run: motion of the animal in N-dimensions
### 2. check.boundary: check whether the boundary of the n-dimensional hyper-cube has been reached. If yes, the animal bounces back from the wall.
### 3. gen.rate: generates the firing rates of N-dimensional PCs in the function of spatial location
### 4. get.rate: calculates the firing rates of N-dimensional PCs in the function of spatial location
### 5. plot.rates: plots the firing rates of N-dimensional PCs in the function of spatial location
### 7. sim.spikes: during theta, neurons fire with he rate defined by the actual position
### during SPW neurons fire according to fantasized movement on the track 15x faster, 5x larger rate

library(viridis)

## a single iteration on the speed
## the speed is an OU process, smoothly changing with a mean of pars.speed$mu.v
## the direction is changed in each timestep by adding a random vector to the current velocity 
## and then renormalising its length
iterate.sv <- function(speed.last, velocity.last, dt, pars.speed){
# dt: time resolution (s)
# pars.speed$mu.v: mean speed (m/s)
# pars.speed$sd.v: sd of the speed 
# pars.speed$tau.v: smoothness of the speed profile (s)
# pars.speed$m.dv: the mean change of velocity - direction
	
	Q <- sqrt(2/pars.speed$tau.v) * pars.speed$sd.v

	speed.new <- speed.last + (pars.speed$mu.v - speed.last) * dt / pars.speed$tau.v + Q * rnorm(1,0,1) * sqrt(dt)
	if (speed.new < pars.speed$v.min)  speed.new <- speed.last + (pars.speed$mu.v - speed.last) * dt / pars.speed$tau.v
    if (speed.new > pars.speed$v.max)  speed.new <- speed.last + (pars.speed$mu.v - speed.last) * dt / pars.speed$tau.v
	
	# change the velocity - random vector of length pars.speed$m.dv
	ndim <- length(velocity.last)
	vv.add <- runif(ndim, -1,1)
	vv.add <- pars.speed$m.dv * vv.add / sqrt(sum(vv.add^2))
	
	velocity.new <- velocity.last + vv.add ## add a random vector
	velocity.new <- velocity.new * speed.new / sqrt(sum(velocity.new^2)) # normalise to speed.new
	sv <- list(speed=speed.new, velocity=velocity.new)
	sv
}


sim.run <- function(Tmax=10, xmax=1, dt=0.01, pars.speed, ndim=2, seed=123, init.v=NULL){
# Tmax <- Tmax; xmax <- 1; dt <- 0.01; ndim <- 2; pars.speed$mu.v <- 0.05; sd.v <- 0.03; 
# pars.speed$tau.v <- 1; pars.speed$m.dv <- 0.02; seed <- 11; pars.speed$v.max <- 0.2; pars.speed$v.min <- 0.02; init.v <- NULL

# # Simulating the movement of the rat on an n-dimensional box
# # the speed is an OU process, smoothly changing with a mean of pars.speed$mu.v
# # the direction is changed in each timestep by adding a random vector to the current velocity 
# # and then renormalising its length
#######################################
# # input:
# # Tmax: duration of the simulations (s)
# # xmax: the size of the n-dimensional square box (m)
# # dt: time resolution (s)
# # ndim: the dimensionality of the random walk
# # pars.speed$mu.v: mean speed (m/s)
# # pars.speed$sd.v: sd of the speed 
# # pars.speed$tau.v: smoothness of the speed profile (s)
# # pars.speed$m.dv: the mean change of velocity - direction
######################################
# # output: a list with elements: 
# # x: position (m)
# # v: velocity (m/s)

	

  	set.seed(seed)
	
	if (pars.speed$mu.v + pars.speed$sd.v > pars.speed$v.max)  stop("mean + sd > v.max")
	if (pars.speed$mu.v - pars.speed$sd.v < pars.speed$v.min)  stop("mean - sd < v.min")
	
	# init the speed
	sv <- list()
	if (is.null(init.v)) sv$speed <- rnorm(1, pars.speed$mu.v, pars.speed$sd.v) else sv$speed <- init.v
	if (sv$speed > pars.speed$v.max)  sv$speed <- pars.speed$mu.v
	if (sv$speed < pars.speed$v.min)  sv$speed <- pars.speed$mu.v
	  
  	#### The hidden state variables
	if (length(xmax) > 1) x.last <- xmax/2 else x.last <- rep(xmax/2, ndim)# we start from the middle
	x.mat <- matrix(NA, ndim, Tmax/dt)
	x.mat[,1] <- x.last
	sv$velocity <- rep(sqrt(sv$speed^2/ndim), ndim) # velocity of leaving
	velocity.mat <- matrix(NA, ndim, Tmax/dt)
	velocity.mat[,1] <- sv$velocity

 	for (i in 1:(Tmax/dt-1)){
		# 1. arriving to the new place
		# print(vv)
	    x.new <- x.last + sv$velocity * dt
		vx <- check.boundary(sv$velocity, x.new, xmax)
		sv$velocity <- vx$v # this does not change the speed only the direction!
		x.new <- vx$x
		x.mat[,i+1] <- x.new
		
		## iterate on the speed and velocity
		sv <- iterate.sv(sv$speed, sv$velocity, dt, pars.speed)
		
		velocity.mat[,i+1] <- sv$velocity
		x.last <- x.new
  	}
  	
  	vx <- list(x=x.mat, v=velocity.mat)
  	vx
}

### check.boundary: check whether the boundary of the n-dimensional hyper-cube has been reached. If yes, the animal bounces back from the wall.

check.boundary <- function(v, x, xmax){
# # input: speed, position, boundaries
# # output: a list with the new speed vector and the new position
# boundary conditions must be fulfilled at each time point!
#  - positions outside the boundary are mapped to the region within the boundary
#  - velocity is negated if the animal moves out
# this works only if boundary is checked at every location!
	if (length(xmax) < length(x)) xmax <- rep(xmax[1], length(x))

	i.turn <- which(((x>xmax) | (x < 0)))
	if (length(i.turn) > 0){
		for (i in i.turn){
			v[i] <- (-1) * v[i]
			if (x[i]>xmax[i]) x[i] <- 2 * xmax[i] - x[i] else x[i] <- (-1) * x[i]
		}
	}
	vx <- list(v=v, x=x)
	vx
}


#################################################################
#################################################################

## generates n-dimensional place cells
gen.rate <- function(N.cells, D.track, volume, scale.length=2, rate0=c(0.1, 0.5), rate.max=c(5, 20), r.place.field=c(0.1, 0.4), alpha=0.57, beta=1/0.14, active.cells = F, graphics = F, seed=317, rand.nfields=T, rand.posfields=T){
	# input:
	# N.cells: number of place cells
	# D.track: dimensionality of the environment
	# volume: the volume of the environment in units of place fields (mean place field diameter, mean(rate.max))
	# scale.length: total volume of the environment scale.length * volume
	# rate0, rate.max: min and max firing rate (Hz)
	# r.place.field: range of radius of place field (m)
	# alpha, beta: parameters of gamma distribution of place field propensity from Rich et al, 2014
	# active.cells: simulate only cells that are active on this environment (at least 1 place fied)
	# graphics: plt the place field
	# rand.nfields: boolean. If True, then number of fields is random. If false, each cell will have one field.
	######################################
	# output:
	# rates - list of N cells with the means and vars of the Gaussian PFs + their baseline firing rate and max rate for each field
	
	unit <- mean(r.place.field) # now in meters
	x.track <- volume ^ (1/D.track) * unit 
	x1.track <- x.track * scale.length
	set.seed(seed)
	
	n.cells.lin <- N.cells ^ (1/D.track)
	
	rates <- list()
	for (i.cell in 1:N.cells){
		if (rand.nfields){
			pf.propensity <- 2 # must be smaller than 1 pfs / m - Rich et al, 2014
			new.cell <- T
			while (new.cell){
				pf.propensity <- rgamma(1, alpha, beta) # place fields / m
				pf.propensity <- pf.propensity * mean(r.place.field) # place fields / unit length
				mu.nfields <- pf.propensity  * volume * scale.length
				if (active.cells) {
					if (mu.nfields > 0.5)  new.cell <- F
				} else {
					new.cell <- F
				}
			}
			n.fields <- 0
			if (active.cells) {
				while (n.fields < 1) n.fields <- rpois(1, mu.nfields)
			} else {
				n.fields <- rpois(1, mu.nfields)
			}
		} else {
			n.fields <- 1
		}
		
		if (n.fields > 0){
			if (rand.posfields) {
				mu <- matrix(runif(n.fields * D.track, 0, x.track), nrow=n.fields)
				mu[,1] <- mu[,1] * scale.length # scale the first dimension according to scale.length
			} else {
				dist.pf <- x.track / (n.cells.lin-3)
				mu <- matrix(c(((i.cell-1) %/% 10) * dist.pf-dist.pf, ((i.cell-1) %% 10) * dist.pf - dist.pf), nrow=n.fields)	
				# mu <- matrix(c(0.1,((i.cell-1) %/% 10) * dist.pf-dist.pf), nrow=1)	
			}
			if (rand.posfields) {
				sdd <- runif(n.fields, r.place.field[1], r.place.field[2])
				sds <- matrix(rep(sdd, D.track), nrow=n.fields)
				# sds <- matrix(runif(n.fields * D.track, r.place.field[1], r.place.field[2]), nrow=n.fields)
				rmax <- runif(n.fields, rate.max[1], rate.max[2])
			} else {
				sds <- matrix(rep(mean(r.place.field), D.track), nrow=n.fields)
				rmax <- rep(mean(rate.max), n.fields)
			}
		} else {
			mu <- matrix(0, nrow=1, ncol=D.track)
			sds <- matrix(0.2, nrow=1, ncol=D.track)
			rmax <- 0
		}
		r0 <- runif(1, rate0[1], rate0[2])
			
		# for (i.place.field in 1:length(mu)) rate <- rate + rr[i.place.field] * exp(-(dist - mu[i.place.field])^2 / (2 * sds[i.place.field]^2))
		rates[[i.cell]] <- list(mu=mu, sds=sds, rmax=rmax, r0=r0)
	}
	attributes(rates) <- list(xmax=c(x.track * scale.length))
	if (D.track > 1) attributes(rates)$xmax <- c(attributes(rates)$xmax, rep(x.track, D.track-1))
	if (graphics) plot.rates(rates)
	rates	
}



########################################
## calculates the firing rate of all cells at a given position
get.rate <- function(xx, rates){
	N <- length(rates)
	ndim <-ncol(rates[[1]]$mu)
	spike.rate <- rep(0, N)
	for (i.cell in 1:N){
		for (i.place.field in 1:length(rates[[i.cell]]$rmax)){
			r.field <- rates[[i.cell]]$rmax[i.place.field]
			for (i.dim in 1:ndim){
				r.field <- r.field * exp(-(xx[i.dim]-rates[[i.cell]]$mu[i.place.field,i.dim])^2/(2*rates[[i.cell]]$sds[i.place.field,i.dim]^2))
			}
			spike.rate[i.cell] <- spike.rate[i.cell] + r.field
		}
		spike.rate[i.cell] <- spike.rate[i.cell] + rates[[i.cell]]$r0
	}
	spike.rate
}

# rates <- gen.rate(9, 2, 64, scale.length=2, graphics=T)
# round(get.rate(c(1,1), rates),2)

# xmax <- 1.6



########################################
## calculates the firing rate of one cell for all positions
##
## in the case of a standard ppc coding (Ma et al., 2006)
## the variance is encoded by the total firing rate of the population: var = 1/r
## the variance to be encoded is between 0.001 and 0.2 
## here I use a scaling factor where at the population level the original rates represent a variance of 0.003 m^2
##
## it is also possible to scale the rate with the sd.
## in this case the var should represent the sd instead and the var.r0 should be the corresponding standard deviation

get.rate.x <- function(xx, rate, variance=NULL, var.r0=0.003){
	L <- ncol(xx)
	ndim <-ncol(rate$mu)
	spike.rate <- rep(rate$r0, L)
	sim.ppc <- F
	if (!is.null(variance)) {
		sim.ppc <- T	
		## variance can be 2D, we can encode only a single uncertainty
		## we take an isotropic Gaussian with the same entropy, that is the geometric mean of the variances
		var.1 <- (apply(variance, 2, prod))^(1/ndim)
	}
	
	for (i.place.field in 1:length(rate$rmax)){
		r.field <- rate$rmax[i.place.field]
		if (sim.ppc) r.field <- r.field * var.r0 / var.1
		for (i.dim in 1:ndim){
			r.field <- r.field * exp(-(xx[i.dim,]-rate$mu[i.place.field,i.dim])^2/(2*rate$sds[i.place.field,i.dim]^2))
		}
		spike.rate <- spike.rate + r.field
	}
	spike.rate
}

# rate <- list(mu=matrix(c(0.5, 0.5), 1), sds=matrix(c(0.15, 0.15),1), rmax=10)
# mu <- rbind(c(0.5, 0.5, 0.5,    0.4, 0.5,    0.4, 0.3, 0.2), c(0.5, 0.5, 0.5,    0.5, 0.4,    0.4, 0.3, 0.2))
# variance <- rbind(c(0.055, 0.1, 0.03, rep(0.055, 5))^2, c(0.055, 0.1, 0.03, rep(0.055, 5))^2)
# get.rate.x(mu, rate)
# get.rate.x(mu, rate, variance)
# get.rate.x(mu, rate, sqrt(variance), var.r0=sqrt(0.003))


########################################
## calculates the firing rate of one cell for all positions
## in the case of a convolutional encoding model (Zemel et al., 1998)
## Product of gaussian densities: Matrix cookbook: 8.1.8
## applied to factorised Gaussians

get.rate.x.int <- function(mu, variance, rate){
	L <- ncol(mu)
	ndim <-ncol(rate$mu)
	spike.rate <- rep(0, L)
	## rate$rmax[i.place.field] is the firing rate if the mean is at the center of the field and the posterior variance is 0.003
	## r_max = c_r \int N(rate$mu, rate$var) N(rate$mu, 0.003) dx
	## c_r = r_max /  \int N(rate$mu, rate$var) N(rate$mu, 0.003) dx 
	## r = c_r * \int N(rate$mu, rate$var) N(mu, variance) dx
	##
	## ONLY if Gaussians are factorised!
	## \int N(rate$mu, rate$var) N(mu, variance) dx = 1/sqrt(det(2*pi*(rate$var + variance))) * \prod_i exp((rate$mu_i - mu_i)^2 / (rate$var_i + variance_i))
	for (i.place.field in 1:length(rate$rmax)){
		r.field <- rate$rmax[i.place.field] * sqrt(prod(2*pi*(rate$sds[i.place.field,]^2 + 0.003))) ## prefactor from c.r
		r.field <- r.field / sqrt(apply(2*pi*(rate$sds[i.place.field,]^2 + variance), 2, prod)) # prefactors from the integral
		for (i.dim in 1:ndim){
			r.field <- r.field * exp(-(mu[i.dim,]-rate$mu[i.place.field,i.dim])^2/(2*(rate$sds[i.place.field,i.dim]^2+variance[i.dim,]))) # exponential factors from the integral
		}
		spike.rate <- spike.rate + r.field
	}
	spike.rate
}



get.rate.conv <- function(mu, variance, rate){
	L <- ncol(mu)
	ndim <-ncol(rate$mu)
	spike.rate <- rep(0, L)
	## rate$rmax[i.place.field] is the firing rate if the mean is at the center of the field and the posterior variance is 0.003
	## r_max = c_r \int N(rate$mu, rate$var) N(rate$mu, 0.003) dx
	## c_r = r_max /  \int N(rate$mu, rate$var) N(rate$mu, 0.003) dx 
	## r = c_r * \int N(rate$mu, rate$var) N(mu, variance) dx
	##
	## ONLY if Gaussians are factorised!
	## \int N(rate$mu, rate$var) N(mu, variance) dx = 1/sqrt(det(2*pi*(rate$var + variance))) * \prod_i exp((rate$mu_i - mu_i)^2 / (rate$var_i + variance_i))
	for (i.place.field in 1:length(rate$rmax)){
		r.field <- rate$rmax[i.place.field] * sqrt(prod(2*pi*(rate$sds[i.place.field,]^2 + 0.003))) ## prefactor from c.r
		r.field <- r.field / sqrt(apply(2*pi*(rate$sds[i.place.field,]^2 + variance), 2, prod)) # prefactors from the integral
		for (i.dim in 1:ndim){
			r.field <- r.field * exp(-(mu[i.dim,]-rate$mu[i.place.field,i.dim])^2/(2*(rate$sds[i.place.field,i.dim]^2+variance[i.dim,]))) # exponential factors from the integral
		}
		spike.rate <- spike.rate + r.field
	}
	spike.rate
}



# # ### Rmax should be the rate when the variance is 0.003
# rate <- list(mu=matrix(c(0.5, 0.5), 1), sds=matrix(c(0.15, 0.15),1), rmax=10)
# mu <- rbind(c(0.5, 0.5, 0.5,    0.4, 0.5,    0.5, 0.5, 0.5), c(0.5, 0.5, 0.5,    0.5, 0.4,    0.3, 0.3, 0.3))
# variance <- rbind(c(0.055, 0.1, 0.03, 0.05, 0.05, 0.1, 0.1, 0.01)^2, c(0.055, 0.1, 0.03, 0.05, 0.05, 0.1, 0.01, 0.1)^2)

# x <- seq(-1, 2, by=0.01)
# xx <- matrix(rep(x, 301), 301)
# yy <- matrix(rep(x, each=301), 301)

# i.post <- 1
# rate.mat <- rate$rmax[1] * exp(-(xx-rate$mu[1,1])^2 / (2*rate$sds[1,1]^2)) * exp(-(yy-rate$mu[1,2])^2 / (2*rate$sds[1,2]^2)) / sqrt(prod(2*pi*rate$sds[1,]^2)) * sqrt(prod(2*pi*(rate$sds[i.place.field,]^2 + 0.003)))
# post.mat <- exp(-(xx-mu[1,i.post])^2 / (2*variance[1,i.post])) * exp(-(yy-mu[2,i.post])^2 / (2*variance[2,i.post])) / sqrt(prod(2*pi*variance[,i.post]))

# sum(rate.mat * post.mat / 100 / 100)

# image(x, x, rate.mat)
# image(x, x, post.mat)
# max(rate.mat)
# max(post.mat)
# sum(post.mat)

########################################
## draws place fields of the cells 
## by default it plots in 1D if 1D data is given and
## plots in 2D otherwise
## unless 3D is explicitely required

plot.rates <- function(rates, plot.3d=F, return.ratemap=F, scale.rate=T){
	N <- length(rates)
	D <- ncol(rates[[1]]$mu)
	xmax <- attributes(rates)$xmax
	if  (D==1){
		x <- seq(0, xmax[1], length=1001)
		rr <- matrix(0, N, 1001)
		for (i.cell in 1:N){
			for (i.place.field in 1:length(rates[[i.cell]]$rmax)){
				rr[i.cell,] <- rr[i.cell,] + rates[[i.cell]]$rmax[i.place.field] * exp(-(x-rates[[i.cell]]$mu[i.place.field])^2/(2*rates[[i.cell]]$sds[i.place.field]^2))
			}
			rr[i.cell,] <- rr[i.cell,] + rates[[i.cell]]$r0
		}
		n.cols <- ceiling(sqrt(N))
		par(mfcol = c(n.cols,n.cols)); par(mar=c(3,3,1,1))
		for (i.cell in 1:N) {
			plot(x, rr[i.cell,], col=rainbow(N)[i.cell], t="l", lty=1, xlab="distance (m)", ylab="firing rate (Hz)", axes=F)
			axis(2, las=2)
		}
		axis(1, padj=-.7)
		mtext("distance (m)", 1, cex=0.7, line=1.5)
	} 
	if ((D > 2) & (plot.3d==T)) {# we are plotting the first two dimentsions
		Lmax <- 20; Lmax1 <- 80 #<- Lmax * xmax[1] / xmax[2]
		rr <- array(0, dim=c(N, Lmax1, Lmax, Lmax))
		for (i.cell in 1:N){
			for (i.place.field in 1:length(rates[[i.cell]]$rmax)){
				ii <- 1
				for (xx in seq(0, xmax[1], length=Lmax1)){
					jj <- 1
					for (yy in seq(0, xmax[2], length=Lmax)){
						kk <- 1
						for (zz in seq(0, xmax[3], length=Lmax)){
							rr[i.cell,ii,jj,kk] <- rr[i.cell,ii,jj,kk] + rates[[i.cell]]$rmax[i.place.field] * exp(-(xx-rates[[i.cell]]$mu[i.place.field,1])^2/(2*rates[[i.cell]]$sds[i.place.field,1]^2)) * exp(-(yy-rates[[i.cell]]$mu[i.place.field,2])^2/(2*rates[[i.cell]]$sds[i.place.field,2]^2)) * exp(-(zz-rates[[i.cell]]$mu[i.place.field,3])^2/(2*rates[[i.cell]]$sds[i.place.field,3]^2))
							kk <- kk + 1
							}
						jj <- jj + 1
					}
					ii <- ii + 1
				}
			}
			rr[i.cell,,,] <- rr[i.cell,,,] + rates[[i.cell]]$r0
			cat(i.cell, " ")
		}
		library(rgl)
		for (i.cell in 1:N) {
			x <- rep(seq(0,xmax[1],length=Lmax1), Lmax^2)
			y <- rep(rep(seq(0,xmax[2], length=Lmax), each=Lmax1), Lmax)
			z <- rep(seq(0,xmax[3], length=Lmax), each=Lmax*Lmax1)
			ps <- as.vector(rr[i.cell,,,])
			ps <- ps - min(ps)
			ps <- ps / mean(ps)
			borders <- quantile(ps, c(0.5, 0.6, 0.7,0.8,0.9, 0.95))
			ps[ps < borders[1]] <- 0
			ps[ps >= borders[6]] <- -6
			ps[ps >= borders[5]] <- -5
			ps[ps >= borders[4]] <- -4
			ps[ps >= borders[3]] <- -3
			ps[ps >= borders[2]] <- -2
			ps[ps >= borders[1]] <- -1
			ps <- ps * (-1)

			x <- x[ps>0]; y <- y[ps>0]; z <- z[ps>0]; ps <- ps[ps>0]
			
			plot3d(x,y,z,pch=21, bg=1, col=rev(rainbow(6, end=0.7))[ps])
			readline(i.cell)


# rr <- array(NA, dim=c(4,3,2))
# for (x in seq(1,4)){
	# for (y in seq(1,3)){
		# for (z in c(1,2)){
			# rr[x,y,z] <- 100*x + 10*y +z			
		# }
	# }	
# }
		
		}
	} else {
		if (D > 1) { # plot in 2D
			Lmax <- 33; Lmax1 <- Lmax * xmax[1] / xmax[2]
			# x <- seq(0, xmax[1], length=Lmax1)
			# y <- seq(0, xmax[2], length=Lmax)
			x <- seq(0+xmax[1]/Lmax1/2, xmax[1]-xmax[1]/Lmax1/2, length=Lmax1)
			y <- seq(0+xmax[2]/Lmax/2, xmax[2]-xmax[2]/Lmax/2, length=Lmax)
			rr <- array(0, dim=c(N, Lmax1,Lmax))
			for (i.cell in 1:N){
				for (i.place.field in 1:length(rates[[i.cell]]$rmax)){
					ii <- 1
					for (xx in seq(0, xmax[1], length=Lmax1)){
						jj <- 1
						for (yy in seq(0, xmax[2], length=Lmax)){
							rmid <- rates[[i.cell]]$rmax[i.place.field] * exp(-(xx-rates[[i.cell]]$mu[i.place.field,1])^2/(2*rates[[i.cell]]$sds[i.place.field,1]^2)) * exp(-(yy-rates[[i.cell]]$mu[i.place.field,2])^2/(2*rates[[i.cell]]$sds[i.place.field,2]^2))
							if (D > 2) {
								for (idim in 3:D){
									rmid <- rmid * exp(-(xmax[idim]/2-rates[[i.cell]]$mu[i.place.field,idim])^2/(2*rates[[i.cell]]$sds[i.place.field,idim]^2))
								}
							}
							rr[i.cell,ii,jj] <- rr[i.cell,ii,jj] + rmid				
							jj <- jj + 1
						}
						ii <- ii + 1
					}
				}
				rr[i.cell,,] <- rr[i.cell,,] + rates[[i.cell]]$r0
				cat(i.cell, " ")
			}
			n.cols <- ceiling(sqrt(N))
			par(mfcol = c(n.cols,n.cols)); par(mar=c(1,1,1,1)/10)
			if (scale.rate){
				brs <- seq(0, ceiling(max(rr)), length=49)
				for (i.cell in 1:N) image(x,y,rr[i.cell,,], br=brs, col=viridis(48), axes=F, xlab='', ylab='')
			} else {
				for (i.cell in 1:N) image(x,y,rr[i.cell,,], col=viridis(48), axes=F, xlab='', ylab='')				
			}
		}
	}
	if (return.ratemap){
		out=list(x=x, y=y, rr=rr)
	} else {
		out=NULL
	}
	out
}

# matplot(t(rates), t="l", col=rainbow(100), lty=1, axes=F, xlab="distance (m)", ylab="firing rate (Hz)"); axis(1); axis(2, las=2)


#################################################################
#################################################################
## simulates the spikes from the location x
## the fantasised location xx (only used when simulating SPW-replay)
## the location dependent firing rate rates
## the state of the system (x or xx is used)
## and the increase in the cell's firing rate during SPWs (gain)
## output:
## NxL matrix with the number of spikes in a given timebin 


sim.spikes.x <- function(x, xx, rates, states, dt, gain.spw=5){
	## x: spatial location on this trial
	## xx: spatial loaction with large speed and extended maze - for fantasizing during sharp waves
	## rates: firing rate as the function of location (Hz)
	## states: maze or sharp vawe?
	## dt: time resolution of x and r (s)
	if (is.vector(x)) {
		x <- matrix(x, nrow=1)
		xx <- matrix(xx, nrow=1)
	}
	
	N.cell <- length(rates)
	L <- ncol(x)
	if (length(states) != L) stop("length of states not equal to the length of positions")
	
	## fantasised trajectories in SPWs
	x[,!!states] <- xx[,!!states]
	gain <- rep(1, L); gain[!!states] <- gain.spw

	## second: generate spike train
	spikes <- matrix(0, N.cell, L)
	for (i.cell in 1:N.cell){
		rate <- rates[[i.cell]]
		spike.rate <- get.rate.x(x, rate)
		spikes[i.cell,] <- rpois(L, spike.rate * dt * gain) 
	}
	spikes
}

# plot(sp[,2], sp[,1], pch=16, cex=0.7, log="")
# plot(sp[-1,2], diff(sp[,2]), t="l", log="y")
# abline(h=0.0004, col=2)

#################################################################
## simulates the spikes from the location x
## same as sim.spikes but hopefully much faster

sim.spikes.x.ppc <- function(mu, variance, rates, dt, conv.enc=T, use_SD=T, sigma0=0.003){
	## mu: represented mean spatial location on this trial
	## var: represented variance of the spatial location on this trial
	## rates: firing rate as the function of location (Hz)
	## dt: time resolution of mu and var (s)
	## conv.enc: boolean; T: convolutional encoding is used. F: standard PPC is used

	if (is.vector(mu)) {
		mu <- matrix(mu, nrow=1)
	}
	if (is.vector(variance)) {
		variance <- matrix(variance, nrow=1)
	}
	
	N.cell <- length(rates)
	L <- ncol(mu)
	
	## generate spike train
	spikes <- matrix(0, N.cell, L)
	for (i.cell in 1:N.cell){
		rate <- rates[[i.cell]]
		if (conv.enc) {
			spike.rate <- get.rate.x.int(mu, variance, rate)
		} else if (use_SD) {
			spike.rate <- get.rate.x(mu, rate, sqrt(variance), var.r0=sqrt(sigma0))
		} else {
			spike.rate <- get.rate.x(mu, rate, variance, var.r0=sigma0)			
		}
		spikes[i.cell,] <- rpois(L, spike.rate * dt) 
	}
	spikes
}


raster2spt <- function(raster){
	t.sp <- NA
	cell.sp <- NA
	N.cells <- nrow(raster)
	for (n.sp in 1:max(raster)){
		for (i.cell in 1:N.cells){
			t.sp <- c(t.sp, rep(which(raster[i.cell,] == n.sp), each=n.sp)) # caution! raster become binary here!
			cell.sp <- c(cell.sp, rep(i.cell, n.sp*length(which(raster[i.cell,] == n.sp))))
		}
	}
	t.sp <- t.sp[-1] / 1000 # convert to seconds
	cell.sp <- cell.sp[-1]
	
	sort.sp <- sort(t.sp, index.ret=T)$ix
	cell.spikes <- cbind(t.sp[sort.sp], cell.sp[sort.sp])
	cell.spikes
}

