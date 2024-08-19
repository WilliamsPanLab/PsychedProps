# Project Overview

This document outlines the steps and methods used in the project. Below is a structured guide for preprocessing, derivations, and analysis.

## 1. Preprocessing

### 1A. fMRI Preprocessing
  fmriprep call - full_fmriprep.sh
  SNR mask derivationi. - TSNR_mask_0_antimask
  SNR mask derivationii. - TSNR_mask_1_ExtractSNR_subjectwise
  SNR mask derivationiii. - TSNR_mask_2_Combine_SNRmaps
  SNR mask derivationiv. - TSNR_mask_3_avSNR-mw_to_label
  SNR mask derivationv. - Downsample
  xcpd call - sbatch_xcpd.sh
  Download psilocybin data - DL_psilo_txts.sh
  Download psilocybin data - DL_psilo_fd.sh
  Download psilocybin data - DL_psilo.sh
### 1B. Ca2+ Preprocessing
  List out sessions: sesh_lister.sh for extracting sessions into text file
  Obtain group-level mask at factor of 2 downsample (NMF) - Group_Mask_and_DS_oneHalf.py
  Obtain group-levcel mask at factor of 6 downsample (OpFl) - Group_Mask_and_DS_oneSixth.py
  Downsample by factor of 2 for NMF - BP_Smooth_oneHalf.py
  Downsample by factor of 6 for Optical flow - DS_Smooth_oneSixth_Drug.py (adapted for each drug with file path and extensions)
### 1C. DMN Derivation
- **1C.I** Human downsampling to fs5
  Downsample script
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
- **1D.V** NMF downsampling
- **1D.VI** DMN validation
  Biol psych
  Yeo7

### 1D. Optical Flow (OpFl) Preprocessing
- **1D.I** Downsampling
  DS humans to fs4
  DS mice to 67x70
- **1D.II** Motion censoring
  Motion Masking
- **1D.III** Running optical flow
  OpFl MDMA
  OpFl Psil
  OpFl Mice

## 2. Optical Flow Derivations

### 2A. Magnitudes
    Extract Magnitudes MDMA
    Extract Magnitudes Psil
    Extract Magnitudes Mice
    Extract Dif *3
### 2B. Bottom-Up Relative Angles
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
