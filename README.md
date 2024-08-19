# Guide to the code behind XXX

This document outlines the steps and methods used in the project. Below is a structured guide for image processing, derivations, and analyses. All image processing was run in a Linux environment using a SLURM cluster for high-compute jobs. In this context, sbatch refers to submitting a job to the SLURM job scheduler. Note that fmriprep and xcpd calls utilize their singularity images, which need to be installed locally.

## 1. Preprocessing

### 1A. fMRI Preprocessing
  This is the [fmriprep](https://fmriprep.org/en/stable/) call used. These derivatives become the inputs of TSNR (temporal signal to noise ratio) scripts as well as xcpd (eXtensible connectivity pipeline) scripts. It is written to run on a SLURM cluster (Simple Linux Utility for Resource Management): - [full_fmriprep](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/full_fmriprep.sh)
  
  After fmriprep has ran, this script is the first in deriving a TSNR mask. It simply masks out the brain from a nifti. This provides an out-of-brain comparison for the magnitude of signal fluctuations: [TSNR_mask_0_antimask](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_0_antimask.sh)
  
  This uses derivatives to calculate SNR for each individual subject: [TSNR_mask_1_ExtractSNR_subjectwise](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_1_ExtractSNR_subjectwise.py)
  
  This script then combines average SNR maps across participants for a global average map. [TSNR_mask_2_Combine_SNRmaps](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_2_Combine_SNRmaps.py)

  Next is a quick downsampling of this mask. We'll need it in fsaverage5 space for regularized non negative matrix factorization. Here's that [script](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DS_surf_TSNR_fs5.sh)

  Next, this script combines the averaged left and right hemisphere TSNR maps. The bottom 15 percent of voxels (lowest TSNR) are masked out for all subsequent analyses. Note that psilocybin TSNR maps indicated that psilocybin study data all met the 15th percentile calculated from the MDMA study, likely due to use of multiecho. The same mask is applied to that dataset for equivalence: [TSNR_mask_3_avSNR-mw_to_label](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/TSNR_mask_3_avSNR-mw_to_label.m)

  Now all we need to do is downsample the mask to appropriate resolution for optical flow, fsaverage4. Here's that [script](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DS_surf_TSNR.sh). Shoutout to [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) (wb_command) for making this vanilla and simple.

  That's it for SNR. Now, let's implement the [xcpd](https://xcp-d.readthedocs.io/en/latest/) sbatching script (useful for head-motion-prone scan protocols, runs on fmriprep outputs): [sbatch_xcpd](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/sbatch_xcpd.sh)

  The WashU crew was kind enough to send their data over already processed. More information on their dataset is available [here](https://www.nature.com/articles/s41586-024-07624-5).For full transparency, the code used to download the images and their associated files are linked below:

  Download psilocybin data, text files: [DL_psilo_txts](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DL_psilo_txts.sh)
  
  Download psilocybin data, framewise displacement files - [DL_psilo_fd](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DL_psilo_fd.sh)
  
  Download psilocybin data: neuroimages - [DL_psilo](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/DL_psilo.sh)

  OK, now we have all our processed data. To make it suitable for NMF and optical flow, we are going to motion-mask it. This goes a step further than normal temporal censoring: we are only interested in sequences of neuroimages uninterrupted by head motion. So rather than simply dropping high-motion frames, we'll only retain sequences of at least 8 TRs that are uninterrupted-by-motion. We'll save out the indices of these sequences for optical flow, as we're only interested in activity movement within clean segments (i.e., if frames 9 and 10 have high motion and we retain frames 1-8 and 11-20, we don't want to try and estimate activity movement between 8 and 11, just within the clean segments.) The script that takes care of this is [MotMask](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/MotMask.m). As will be the case for many scripts below, an [equivalent script](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/MotMask_psil.m) exists for the psilocybin study. In general, the only differences are the input filepath and output file path. For this script in particular, some additional lines of code exist. As we get further in the pipeline, the file derivatives will be more equivalent, requiring fewer differences between MDMA study and Psilocybin Study scripts.

  Nice. Done with fMRI preprocessing!
  
### 1B. Ca2+ Preprocessing
  The first thing to do is to list out sessions to process. This will enable us to loop over sessions to-be-processed within python: [sesh_lister](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/sesh_lister.sh). Note that several different sesh listers exist for different drugs and pre vs. post conditions, but all use the same framework. More information on this dataset is available [here](https://www.nature.com/articles/s41586-020-2731-9).

  Next, we'll obtain a group-level mask so that no mouse has pixels included that other mice don't. We'll combine this with downsampling to saveout a mask that is applicable to processing in downsampled space. As for human data, there are two different downsampled resolutions, with the greater resolution being provided to NMF and the lesser resolution being provided to optical flow for computational tractability.

  Obtain group-level mask at factor of 2 downsample (NMF) - [Group_Mask_and_DS_oneHalf](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Group_Mask_and_DS_oneHalf.py)

  Obtain group-levcel mask at factor of 6 downsample (OpFl) - [Group_Mask_and_DS_oneSixth](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/Group_Mask_and_DS_oneSixth.py)

  Just for cleanliness, these scripts are ran independently of downsampling the data for use in NMF and optical flow (even though internal calculations are matched). Here are those scripts. Note they also include gaussian smoothing [for SNR](https://support.brainvoyager.com/brainvoyager/functional-analysis-preparation/29-pre-processing/86-spatial-smoothing) purposes.

  Downsample by factor of 2 for NMF - [BP_Smooth_oneHalf](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/BP_Smooth_oneHalf.py)

  Downsample by factor of 6 for Optical flow - [DS_Smooth_oneSixth_Drug](https://github.com/WilliamsPanLab/PsychedProps/blob/master/scripts/mice/BP_Smooth_oneHalf.py) (adapted for each drug with file path and extensions)

   By the end of this section, you might notice that the mouse data is generally a little less cumbersome to process. That pattern will continue. 

### 1C. DMN Derivation
Before evaluating how DMN function is altered, we have to derive our DMN definition. This is done with regularized non negative matrix factorization. See [this link](https://www.sciencedirect.com/science/article/pii/S1053811917303944?via%3Dihub) for the original paper. 

- **1C.I** Human downsampling to fs5
  First, the ciftis are downsampled to fsaverage5 to play nicely with NMF. That downsample script is available [here]()
- **1C.II** Human NMF: extra scans and NMF parameters
  NMF script 1
  NMF script 2
  NMF script 3
- **1C.III** Mouse NMF: algorithm alterations, extra scans, NMF parameters
  NMF script 1
  NMF script 2
  NMF script 3
- **1C.IV** NMF thresholding
  Script picking out DMN humans
  Script picking out DMN mice
- **1C.V** NMF downsampling

### 1D. DMN validation
  Biol psych
  Yeo7


## 2. Optical Flow Derivations
### 2A. Running optical flow
  OpFl MDMA
  OpFl Psil
  OpFl Mice
  
### 2B. Magnitudes
  Extract Magnitudes MDMA
  Extract Magnitudes Psil
  Extract Magnitudes Mice
  Extract Dif *3
  
### 2C. Bottom-Up Relative Angles
  Extract Relative angles MDMA
  Extract Relative angles Psil
  Extract Relative angles mice
  Extract Dif *3
    
### 2C. Contrast Maps
- **2C.I** Humans
  Aggregate AvgMags
  Aggregate AvgBups
- **2C.II** Mice
  Aggregate AvgMags
  Aggregate AvgBups

## 3. DMN Measurement Derivations

### 3A. Functional Connectivity (FC)
- **3A.I** Humans
  Extract_DMNSeg
  Extract_DMNSeg_psil
  Extract dif *2
- **3A.II** Mice
  Extract_DMNSeg
  Extract dif
### 3B. Autocorrelation
- **3B.I** Humans
  Extract_DMNTA
  Extract_DMNSeg_psil
  Extract dif *2
- **3B.II** Mice
  Extract_TAutoCor
  Extract_TAutoCor_psil
  Extract dif
## 4. Main Effects
### 4A. DMN integration
- **4A.I** MDMA
  .rmd
- **4A.II** Psilocybin
  .rmd (note this will also prep saveouts for unified DMN analyses)
- **4A.III** LSD
  .rmd
### 4B. DMN Autocor
- **4B.I** MDMA
  .rmd
- **4B.II** Psilocybin
  .rmd (note this will also prep saveouts for unified DMN analyses)
- **4B.III** LSD
  .rmd
### 4C. Magnitudes
- **4C.I** MDMA
  .rmd
- **4C.II** Psilocybin
  .rmd (note this will also prep saveouts for unified DMN analyses)
- **4C.III** LSD
  .rmd

### 4D. Bottom-Up Analysis
- **4D.I** MDMA
  .rmd
- **4D.II** Psilocybin
  .rmd
  .rmd of lasting effects
- **4D.III** LSD
  .rmd

## 5. Bootstraps and AUC curves
.rmd MDMA
.rmd psil
.rmd mice

## 6. Self-Report: MDMA
.rmd MDMA
## 7. Brain Visualizations
  Viz DMN grad humans
  Viz DMN grad mice
  Viz Opfl vectors humans
  Viz Opfl vectors mice
  Connectome Workbench visualizations
