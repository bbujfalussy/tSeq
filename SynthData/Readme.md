Computer code associated with our manuscript *Sampling motion trajectories during hippocampal theta sequences*, eLife, 2022.

author: Balazs B Ujfalussy (ujfalussy.balazs@koki.hu)

All scripts use R (version 4.1.2 (2021-11-01) -- "Bird Hippie")


This folder contains code that
- Generates synthetic place cell activity data with embedded theta sequences encoding hypothetical motion trajectories using one of four possible coding schemes.
- Analyse synthetic place cell activity

Our simulations involve generating random trajectories and random spikes. As random seeds are not set, rerunning the simulations will result in slight differences in the individual datapoints.


Code structure:


## I. Generation of synthetic data:
	main_sim.R: 
		Generate synthetic theta-trajectories and 
		encode them in simulated place cell activity using probabilistic coding schemes (used for Figures 3-6)

		Output: ./motiondata_*.RData files 
			./spikes_*.RData files 


	sim_thetaseq_kalman_1D.R: generating theta sequences in 1D (Figure 4–Figure supplement 3.)

	sim_thetaseq_kalman_cycling.R: generating theta sequences with specified autocorrelation structure (as in Fig. 7 a-f)

## II. Analysis of synthetic data
	main_anal.R: The main function to analyse the data
		controls preprocessing and performs various analysis on data from different schemes.
		For more information, see comments in the file 

		Output: Generated figures are saved in subdirectories corresponding to Figures of the paper in the directory SimFigures/
			The subfolders are created automatically.
			analysis results are often saved into RData files in the main folder or in the subfolders
			This can speed up the analysis when repeated multiple times.

## Scripts loaded by the main control script

	preproc_syntetic_data.R: Preprocessing synthetic data - calculating ratemaps and decoding position in three temporal bins 
		using equal bin duration or equal spike counts (if evenspikes==TRUE) within theta cycles. 

	calc_plot_8equalDur.R: calculates the (forward and lateral) bias of the decoder and the firing rate as a function of theta phase
		- using standard ratemaps
		- 8 temporal windows of 120deg width each

	decode_DDC.R:  DDC decoding using standard ratemaps

	SamplingIndex2.R: calculates cycle-to-cycle variability and the sampling index

## Other files, implementing independent analysis pieces

	Combine_DDC.R: combines the analysis performed on the various coding schemes 
		and creates the plots of the results 
		replicating Fig 5b,c,d

	anal_thetaseq_kalman_CCP.R: Calculate CCP - the length of constant cycling periods for synthetic data
		output: SimFigs/10_CorrelatedSamples/

## Other code, defining functions not to be run independently:

	SimRunND.R: Functions used to generate movement trajectory in N-Dimensions and position-dependent neuronal activity

	gen_inf_kalman.R: functions to generate data and perform inference under the Kalman filter model detailed in Figure 3–Figure supplement 1a. 
		In short, the agent is trying to follow an extarnally set goal path, using noisy motor commands and its noisy position estimate

	utils.R: minor utils helpful for plotting or data transformations

	decoding_curves.R: This is a script that calculates the mean decoding error as a function of the position shifted in space or in time. 
		Its output is similar to Fig 2d, Fig 3f 

		It is used during preprocessing and also for calculating TEE when testing sampling versus mean schemes.



