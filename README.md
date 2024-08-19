# Project Overview

This document outlines the steps and methods used in the project. Below is a structured guide for preprocessing, derivations, and analysis.

## 1. Preprocessing

### 1A. fMRI Preprocessing

### 1B. Ca2+ Preprocessing

### 1C. Optical Flow (OpFl) Preprocessing
- **1C.I** Downsampling
- **1C.II** Motion censoring
- **1C.III** Time-series reconstruction in "good frame" segments
- **1C.IV** Running optical flow

### 1D. DMN Derivation
- **1D.I** Human downsampling to fs5
- **1D.II** Human NMF: extra scans and NMF parameters
- **1D.III** Mouse NMF: algorithm alterations, extra scans, NMF parameters
- **1D.IV** NMF thresholding
- **1D.V** NMF downsampling

## 2. Optical Flow Derivations

### 2A. Magnitudes

### 2B. Bottom-Up Relative Angles

### 2C. Contrast Maps
- **2C.I** Humans
- **2C.II** Mice

## 3. DMN Measurement Derivations

### 3A. Functional Connectivity (FC)
- **3A.I** Humans
- **3A.II** Mice

### 3B. Autocorrelation
- **3B.I** Humans
- **3B.II** Mice

## 4. Main Effects

### 4A. Magnitudes
- **4A.I** MDMA
- **4A.II** Psilocybin
- **4A.III** LSD

### 4B. Bottom-Up Analysis
- **4B.I** MDMA
- **4B.II** Psilocybin
- **4B.III** LSD

## 5. Bootstraps

## 6. AUC Curves

## 7. Self-Report: MDMA
## 8. Brain Visualizations
