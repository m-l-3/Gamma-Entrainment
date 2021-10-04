# Overview
In this repository we provide all the implemented codes and materials used for investigating several temporal/spatial synchronization features as well as cross-frequency coupling effects in the brain oscillatory activity recorded through electroencephalography (EEG). The EEG data are included in the *Auditory data*, *Rest data*, and *Multisensory data* folders. 

To see data specifications please see [Auditory Gamma Entrainment](https://openneuro.org/datasets/ds003800/versions/1.0.0) and [Multisensory Gamma Entrainment](https://openneuro.org/datasets/ds003805/versions/1.0.0) for more information on how the data are acquired and what each dataset contains.

The attached MATLAB scripts are the codes that were used to generate the plots and produce the analytical results discussed in the article ***Gamma Entrainment Improves Synchronization Deficits Caused by Dementia*** authored by M. Lahijanian, H. Aghajan, Z. Vahabi, and A. Afzal. The preprint version of the article is accessible through https://doi.org/10.1101/2021.09.30.462389.

# Requirements
For running the scripts, you only need MATLAB. In addition, for producing the connectivity graphs as in Fig. 4b of the referenced article, we recommend using [BrainNet Viewer](https://www.nitrc.org/projects/bnv/).

# Instructions
There are two main scripts: `pre.m` and `main.m` in the [*Code*](/Code) folder. Initial parameters are set at the beginning of each script. First, to define epochs in the data and perform preprocessing steps, run `pre.m` for different window lengths (change `W` to the desired window length). This will prepare the data for the following analyses and the result should be saved in the *Partitioned data* folder (you need to make *Partitioned data* folder in the main repository). Then to generate the epoched data try running `pre.m` for window lengths of 1, 20, and 40 seconds since these values are used later in the analyses. Pay attention to the commented notes as you need to comment/uncomment some lines of the script for `W = 40`.

Then, for generating the results go through `main.m` and you can see the script is pretty segmented to different parts and each part is related to each result discussed in the paper. Try to run each segment one after another. At the beginning of each segment, the epoched data is loaded and some parameters are set to desired values. Then, the segment will end by generating figures as depicted in the paper. To regenerate all the sub-figures you may need to change some parameters and comment/uncomment some lines of the script. Pay attention to the commented notes. 

Note that running cross frequency coupling analyses in `Theta-Gamma Coupling` and `Audio + Visual Entrainment` segments may take a considerable amount of time and demands extensive CPU and memory capacity to generate the coupling features from the epoched data. To avoid running this part, you can used the pre-calculated values of these features which are enclosed in the [*Code*](/Code) folder. This allows you to skip running the calculation parts and directly go to the visualization segments.

# Citation
In the case of using the data in this repository, pleas cite https://doi.org/10.1101/2021.09.30.462389.

# References
Some pre-existed scripts and functions have been used in different parts of `main.m` which are referenced below:
- EEGLab2020.0. https://sccn.ucsd.edu/eeglab/index.php
-	Charles (2021). cbrewer : colorbrewer schemes for Matlab 
(https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab), MATLAB Central File Exchange. Retrieved September 7, 2021.
-	Víctor Martínez-Cagigal (2021). Topographic EEG/MEG plot
(https://www.mathworks.com/matlabcentral/fileexchange/72729-topographic-eeg-meg-plot), MATLAB Central File Exchange. Retrieved September 7, 2021.
-	Felix Zoergiebel (2021). ploterr 
(https://www.mathworks.com/matlabcentral/fileexchange/22216-ploterr), MATLAB Central File Exchange. Retrieved September 7, 2021.
-	Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project
https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847
-	Simon Musall (2021). stdshade 
(https://www.mathworks.com/matlabcentral/fileexchange/29534-stdshade), MATLAB Central File Exchange. Retrieved September 7, 2021.
- BrainNet Viewer (version 1.7). Xia, M., Wang, J. & He, Y. BrainNet Viewer: A Network Visualization Tool for Human Brain Connectomics. PLoS One 8, 68910 (2013).
-	https://github.com/anne-urai/Tools/blob/master/plotting/mysigstar.m
-	https://github.com/muntam/TF-PAC
