# Project Overview

This document outlines the steps and methods used in the project. Below is a structured guide for preprocessing, derivations, and analysis.

## 1. Preprocessing

### 1A. fMRI Preprocessing
fmriprep call
SNR mask derivation

### 1B. Ca2+ Preprocessing
BP and mask script

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
