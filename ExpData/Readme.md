Computer code associated with our manuscript *Sampling motion trajectories during hippocampal theta sequences*, eLife, 2022.

author: Balazs B Ujfalussy (ujfalussy.balazs@koki.hu)

All scripts use R (version 4.1.2 (2021-11-01) -- "Bird Hippie")

This folder contains code that perform analyses of the data from Pfeiffer and Foster, Hippocampal place-cell sequences depict future paths to remembered goals. Nature. 2013 May; 497(7447):74–79. 

Currently we don!t have the permission to publish the data. To obtain the data, please contect the corresponding author of the original paper.

Our analysis involves distributing spikes into different theta phases (e.g., early, mid and late) with equal spike counts. Cycles with spike counts not multiple of 3 are divided randonmly. As random seeds are not set, rerunning the simulations will result in slight differences in the individual datapoints.


## Main analysis control

main_anal.R: The main function to analyse the data
	controls preprocessing and performs various analysis on all sessions.
	For more information, see comments in the file 

	Output: DataFigures/
		Generated figures are saved in subdirectories corresponding to Figures of the paper in the directory DataFigures/
			The subfolders are created automatically.
			analysis results are often saved into RData files in the main folder or in the subfolders
			This can speed up the analysis when repeated multiple times.

## Scripts loaded by the main control script

	preprocess_ephys_data.R
		process motion data, load LFP, detect LFP artefacts and calculate ratemaps
		also reproduces Fig 3 Figure supplement 2. The motion profile of the animal in the experimental sessions. 
		and Figure 3–Figure supplement 3. Place cell firing in the experimental data. 

	calc_plot_seqs.R
		decode theta sequences - 20ms windows with 5ms delay
		reproduces Figure 2c,f,g

	calc_plot_8equalDur.R
		in each 120-deg window shifted by 45-deg calculate the firing rate and the decoding bias
		reprodices Figure 4 e-h; Figure 4–Figure supplement 2

	MAP_decode.R
		select segments from the early, mid and late phase of the theta with equal number of spikes
		Figure 4–Figure supplement 2: decoding error as a function of spatially or temporally shifted real position (see also Fig. 2D)
		Fig 2E: place field size estimated from early and late spikes 


	DDC_decode.R
		Calculates the DDC decoded posterior - as in Fig 5e.

	calc_sampling_index2.R
		calculating CCV and TEE - as in Fig 6h-j and Figure 6–Figure supplement 1g,h,i

	calc_sample_acf.R
		calculating the cycling length and the cycling index (Fig. 7g-i).

## Other files, implementing independent analysis pieces

	ViewData.R: Loads the spike data

	ViewLFPData.R:  Loads the LFP data

	TestMultipleFrames.R: Calculates the ratemaps from home and away trials separately. 
		Also calculates shuffle controls by reconstructing ratemaps from random partiotionning of the data.

		Reproduces Figure6–Figure supplement 2a-d

		Output files: DataFigures/Fig6_FS2/

	motion_sampling.R:
		Generates synthetic data based on potential motion trajectories sampled from aligned real trajectories of the animal - as in Figure 7–Figure supplement 2.
		Saves data into DataFigures/Fig7_FS2/spikes_.RData
		The same analysis pipeline can be run on this data as on the synthetic data (../SynthData)


## Scripts and useful functions
	
	bandFilt.R: quick and easy filter based on fast Fourier transformation
		used to filter LFP in the theta band and detect artefacts

	decoding_curves_exp.R: decoding error as a function of spatially or temporally shifted real position (see also Fig. 2D)


	PlaceCellFunctions.R: 
		Functions related to analysis of place cell activity

	SpeedFunctions.R
		Functions related to analysing motion trajectories





