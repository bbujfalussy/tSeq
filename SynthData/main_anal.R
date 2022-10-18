## Decoding synthetic hippocampal theta sequences. 

##############################################
##############################################
# # 1. MAP or mean scheme: no uncertainty is represented in the population activity (PA)
# # 		- PA represents the the mean
# # 		- the error between two subsequent trajectories is SMALLER than the error between the true and a single trajectory
# # 		- trajectories are necoded with CONSTANT precision
# #
# # 2. sampling scheme
# #		- PA samples from the posterior
# # 		- the error between two subsequent trajectories is SIMILAR to the error between the true and one of the trajectories
# # 		- trajectories are necoded with CONSTANT precision
# #
# # # 3. product scheme (we call it PPC for historical reasons)
# # 		- the PA encodes the full posterior at any given time
# # 		- the error between two subsequent trajectories is SMALLER than the error between the true and a single trajectory
# # 		- trajectories are necoded with DECREASING precision
# #		- variance is encoded in the firing rate
# #
# # # 4. DDC scheme
# # 		- the PA encodes the full posterior at any given time
# # 		- the error between two subsequent trajectories is SMALLER than the error between the true and a single trajectory
# # 		- trajectories are necoded with DECREASING precision
# # 		- variance is encoded in the population activity


## 1. load the required packages
graphics <- F
set.seed(33)
library(Matrix)
library(viridis)
library(shape)
require(gplots)
require(colormap)

source('./gen_inf_kalman.R', chdir = TRUE)
source('./utils.R', chdir = TRUE)
source('../ExpData/PlaceCellFunctions.R')
source('../ExpData/SpeedFunctions.R', chdir = TRUE)


Tmax <- 300 # 300 or 1800
N.cells <- 200
codes <- c('MAP', 'sampling', 'DDC', 'PPC_var')
irregs <- c('', 'irreg')
jitters <- c(0, 5, 10, 20, 30, 40) # ms, uniform random noise on theta start and end times
samplers <- c('', '_same', '_similar', '_uniform', '_anti-quarter', '_anti-half', '_strong-anti')

i.code <- 1 # 1-4
i.irreg <- 1 # 1-2
i.jitter <- 1 # 1-6, we only used it with sampling and MAP code 
i.sampler <- 1

spPth <- 32
estimrates <- TRUE
evenspikes <- TRUE
if (evenspikes) evensp <- '_evenspikes' else evensp <- ''
savefile <- F

for (i.irreg in 1:2){
# i.irreg <- 2
	for (i.code in c(3,4,5)){
		jitmax <- 1
		if ((i.irreg == 2)&(i.code < 3)) jitmax <- 6
		for (i.jitter in 1:jitmax){
			# for (i.sampler in 2:7){
			for (i.sampler in c(1)){
				# i.jitter <- 3
				code <- codes[i.code]
				reg <- irregs[i.irreg]
				jitter <- jitters[i.jitter]
				sampler <- samplers[i.sampler]
				
				###########################################
				## preprocessing:
				## - ratemaps
				## - MAP decoding
				## - decoding shifts
				## - maps based on spikes in the early phase of theta
				## - sharpmaps
				###########################################
				# if (i.code == 5) calc_sharpmaps <- T
				source('preproc_syntetic_data.R')	
	
				###########################################
				## calculate the (forward and lateral) bias of the decoder and the firing rate as a function of theta phase
				## - using standard ratemaps
				## - 8 temporal windows of 120deg width each
				
				graphics <- T
				source('calc_plot_8equalDur.R')
			
				###########################################
				## DDC decoding
				## - using standard ratemaps
				###########################################
				
				source('decode_DDC.R')
								
				###########################################
				## sampling index
				###########################################

				source('SamplingIndex2.R')
				cat('\n \n', i.code, i.jitter, '\n')
			}
		}
	}
}
	
