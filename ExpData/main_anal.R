########################################################################
## combine LFP, position and spike data to estimate the decoding bias at different phases of the theta cycle - 
## find the average starting phase of the theta sequences

library(R.matlab)
library(colormap)
library(ellipse)
library(gplots)
library(seewave)

source('./SpeedFunctions.R')
source('./PlaceCellFunctions.R')
source('./bandFilt.R', chdir = TRUE)
source('../SynthData/utils.R')


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
i_day <- 2

## theta phase that corresponds to the boundary of theta cycles
## these numbers were calculated from the cosine fit to the forward bias - as in Decode_theta_phase - theta_bias.pdf
theta_peaks <- matrix(NA, 2, 4)
theta_peaks <- matrix(c(2.32, 2.27, 2.41, 1.11, 1.9, 1.98, 1.64, 1.68), 2, 4)
graphics <- T
savefile <- 0
dir.create('./DataFigures/', showW=F)


for (i_rat in c(1,2,3,4)){
	start_day <- 1
	end_day <- 2
	# if (i_rat == 1) start_day <- 2
	# if (i_rat == 1) end_day <- 1
	# if (i_rat == 3) start_day <- 2
	for (i_day in start_day:end_day){
		##################################################################################
		## process motion data, load LFP, detect artefacts and calculate ratemaps
		## also reproduces Fig 3 Figure supplement 2. The motion profile of the animal in the experimental sessions. 
		## and Figure 3–Figure supplement 3. Place cell firing in the experimental data. 
		##################################################################################
		graphics <- F
		source('preprocess_ephys_data.R')		
		theta_peak <- theta_peaks[i_day, i_rat]
		theta_start <- theta_peak + 3*pi/4
		
		##################################################################################
		# ## decode theta sequences - 20ms windows with 5ms delay - Figure 2c,f,g
		# ########################################################################

		source('calc_plot_seqs.R')

		##################################################################################
		## in each 120-deg window shifted by 45-deg calculate the firing rate and the decoding bias
		##  - Figure 4 e-h; Figure 4–Figure supplement 2
		##################################################################################

		source('calc_plot_8equalDur.R')	

		## identify the start of the theta sequences 
		## these numbers are from the peak of phase from the Decode_theta_phase - theia_bias.pdf
		theta_peak <- newpars$par[3]
		theta_peaks[i_day, i_rat] <- theta_peak

		# ##################################################################################
		# ## select segments from the early, mid and late phase of the theta with equal number of spikes
		# ########################################################################
		# Fig. 2D and Figure 4–Figure supplement 2: decoding error as a function of spatially or temporally shifted real position 
		# 			DataFigures/05_decoding_curves/
		# Fig 2E: place field size estimated from early and late spikes 
		# 			DataFigures/07PFsize
		source('MAP_decode.R')

		# ##################################################################################
		# ## decode using DDC...
		# ########################################################################
		# Calculates the DDC decoded posterior - as in Fig 5e.
		
		source('DDC_decode.R')

		# ##################################################################################
		# ## calculate the sampling index
		# ########################################################################
		# calculating CCV and TEE - as in Fig 6h-j and Figure 6–Figure supplement 1g,h,i

		for (include_trials in c('home', 'away', 'all')){
			source('calc_sampling_index2.R')
		}
		
		## calculating the cycling length and the cycling index (Fig. 7g-i).
		source('calc_sample_acf.R')
		cat('rat', i_rat, ' day', i_day, ' finished', sep='')
	}
}
