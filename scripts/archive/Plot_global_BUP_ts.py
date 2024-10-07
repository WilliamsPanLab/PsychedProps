# import needed stuff
import scipy.io as sio
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
matplotlib.use('Agg')
import numpy as np
import sys
import os
import pandas as pd
from scipy.signal import find_peaks
# sliding window for plots
def sliding_window_average(data, window_size=1):
    window = np.ones(window_size) / window_size
    return np.convolve(data, window, mode='valid')

# load in subject session dose correspondence
subSeshDose=pd.read_table('/home/users/apines/subjSeshDoseCorresp.csv',header=None,delimiter=' ')

# get number of subjects
numSubjs=subSeshDose.shape[0]

# get placebo sessions
placSeshs=subSeshDose[2]
# 80
mg80Seshs=subSeshDose[3]
# 120
mg120Seshs=subSeshDose[4]

# Loop over each subject
for idx in [0,1,2,4,5,6,8,10,11,12,13,14,15,16]:
for idx in [15,16]:
    # for resting state 1
    task = 'rs1'
    # Extract subject and session info
    subj = subSeshDose.iloc[idx, 0]
    subj = f'sub-MDMA{int(subj):03d}'
    placSesh = placSeshs[idx] 
    m80Sesh = mg80Seshs[idx]
    m120Sesh = mg120Seshs[idx]
    # initialize blank time series to prevent plotting issues later
    plac_ts=0
    m80_ts=0
    m120_ts=0
    # Directory path
    childfp = f'/scratch/users/apines/data/mdma/{subj}'
    filepath_l = f'{childfp}/{placSesh}/{subj}_{placSesh}_{task}_k1_Prop_TS_dmn_L.csv'
    filepath_r = f'{childfp}/{placSesh}/{subj}_{placSesh}_{task}_k1_Prop_TS_dmn_R.csv'
    if os.path.exists(filepath_l) and os.path.exists(filepath_r):
        aTS_l = pd.read_csv(filepath_l).values
        aTS_r = pd.read_csv(filepath_r).values
        combined_ts = np.vstack((aTS_l,aTS_r))
        plac_ts = np.mean(combined_ts,0)
        plac_ts = sliding_window_average(plac_ts)
    # 80 mg
    filepath_l = f'{childfp}/{m80Sesh}/{subj}_{m80Sesh}_{task}_k1_Prop_TS_dmn_L.csv'
    filepath_r = f'{childfp}/{m80Sesh}/{subj}_{m80Sesh}_{task}_k1_Prop_TS_dmn_R.csv'
    if os.path.exists(filepath_l) and os.path.exists(filepath_r):
        aTS_l = pd.read_csv(filepath_l).values
        aTS_r = pd.read_csv(filepath_r).values
        combined_ts = np.vstack((aTS_l,aTS_r))
        m80_ts = np.mean(combined_ts,0)
        m80_ts = sliding_window_average(m80_ts)
    # 120 mg
    filepath_l = f'{childfp}/{m120Sesh}/{subj}_{m120Sesh}_{task}_k1_Prop_TS_dmn_L.csv'
    filepath_r = f'{childfp}/{m120Sesh}/{subj}_{m120Sesh}_{task}_k1_Prop_TS_dmn_R.csv'
    if os.path.exists(filepath_l) and os.path.exists(filepath_r):
        aTS_l = pd.read_csv(filepath_l).values
        aTS_r = pd.read_csv(filepath_r).values
        combined_ts = np.vstack((aTS_l,aTS_r))
        m120_ts = np.mean(combined_ts,0)
        m120_ts = sliding_window_average(m120_ts)
    # get combined time series to caclulate mean and SD
    combinedTS=np.hstack((plac_ts,m80_ts,m120_ts))
    combinedMean=np.mean(combinedTS) 
    print(combinedMean)
    combinedSD=np.std(combinedTS)
    print(combinedSD)
    placPeaks, _ = find_peaks(plac_ts,height=combinedMean+(3*combinedSD))
    m80Peaks, _ = find_peaks(m80_ts,height=combinedMean+(3*combinedSD))
    m120Peaks, _ = find_peaks(m120_ts,height=combinedMean+(3*combinedSD))
    plt.figure(figsize=(10, 6))
    plt.plot(plac_ts, color='black', alpha=0.5, label='Placebo')
    plt.plot(placPeaks, plac_ts[placPeaks], 'ko')
    plt.plot(m80_ts, color='red', alpha=0.5, label='80 mg')
    plt.plot(m80Peaks, m80_ts[m80Peaks], 'ro')
    plt.plot(m120_ts, color='blue', alpha=0.5, label='120 mg')
    plt.plot(m120Peaks, m120_ts[m120Peaks], 'bo')
    plt.title(f'Subject {subj} - {task}')
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.savefig(f'{childfp}/{task}_Subject_{subj}.png')
    plt.close()

