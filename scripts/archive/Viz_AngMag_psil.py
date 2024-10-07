import scipy.io as sio
import nibabel as nb
#import hcp_utils as hcp
import numpy as np
import nilearn.plotting as plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import copy
import sys
import os
from PIL import Image

# import arguments
subj = sys.argv[1]
sesh = sys.argv[2]
task = sys.argv[3]
# import magnitude/angle time series
MagAngL_fp='/scratch/users/apines/data/psil/' + subj + '/' + sesh + '/' + subj + '_' + sesh + '_task-' + task + '_p2mm_masked_Bold_and_Angles_L.mat'
MagAngR_fp='/scratch/users/apines/data/psil/' + subj + '/' + sesh + '/' + subj + '_' + sesh + '_task-' + task + '_p2mm_masked_Bold_and_Angles_R.mat'
MagAng_L=sio.loadmat(MagAngL_fp)['FullMatrix_L']
MagAng_R=sio.loadmat(MagAngR_fp)['FullMatrix_R']
# isolate magnitude and angles
Mag_L=MagAng_L[:,:,0]
Ang_L=MagAng_L[:,:,1]
Mag_R=MagAng_R[:,:,0]
Ang_R=MagAng_R[:,:,1]
Ang_L_rad = np.radians(Ang_L)
Ang_R_rad = np.radians(Ang_R)

# Define the number of angle bins and magnitude bins
num_angle_bins = 10  # 0-180 degrees divided into 18 bins (10 degrees each)
num_magnitude_bins = 10

# Create angle bins and magnitude bins
angle_bins = np.linspace(0, np.pi, num_angle_bins+1)  # Adjust for 0-180 degrees
# Define percentiles for magnitude bins (adjust as needed)
percentiles = np.linspace(0, 100, num_magnitude_bins+1)  # 0-100% divided into num_magnitude_bins percentiles
magnitude_percentiles_L = np.percentile(Mag_L, percentiles)
magnitude_percentiles_R = np.percentile(Mag_R, percentiles)

# Bin the data using global percentiles for each angle bin: right
binned_percentile_magnitudes_L = np.zeros((num_angle_bins, num_magnitude_bins))
for i in range(num_angle_bins):
    angle_min = angle_bins[i]
    angle_max = angle_bins[i + 1]
    mask = (Ang_L_rad >= angle_min) & (Ang_L_rad < angle_max)
    # Calculate percentile bins for the magnitude data using global percentiles
    bin_edges = magnitude_percentiles_L
    # Create histograms using the calculated bin edges
    binned_percentile_magnitudes_L[i, :], _ = np.histogram(Mag_L[mask], bins=bin_edges)

# Bin the data using global percentiles for each angle bin: right
binned_percentile_magnitudes_R = np.zeros((num_angle_bins, num_magnitude_bins))
for i in range(num_angle_bins):
    angle_min = angle_bins[i]
    angle_max = angle_bins[i + 1]
    mask = (Ang_R_rad >= angle_min) & (Ang_R_rad < angle_max)
    # Calculate percentile bins for the magnitude data using global percentiles
    bin_edges = magnitude_percentiles_R
    # Create histograms using the calculated bin edges
    binned_percentile_magnitudes_R[i, :], _ = np.histogram(Mag_R[mask], bins=bin_edges)

# combine percentile bins
binned_percentile_magnitudes=binned_percentile_magnitudes_L+binned_percentile_magnitudes_R;
# save out for group-level merging
output_csv_path='/scratch/users/apines/data/psil/' + subj + '/' + sesh + '/' + subj + '_' + sesh + '_task-' + task + '_DMN_2dhist.csv'
np.savetxt(output_csv_path, binned_percentile_magnitudes, delimiter=',')

# Create a polar contour plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

# define extent (note radial axis is in radians, but labeled as degrees for interp.)
extent = [0, np.pi, 0, 100]
# Create polar plot using contourf
contour = ax.contourf(binned_percentile_magnitudes.T,extent=extent,cmap='inferno')

# Add a colorbar
cbar = plt.colorbar(contour, ax=ax)
cbar.set_label('Percentile Magnitude')

# Customize the plot if needed (e.g., labels, title, etc.)
ax.set_thetamin(0)  # Set the starting angle
ax.set_thetamax(180)  # Set the ending angle
ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
ax.set_theta_zero_location("N")  # Set the zero angle to the north
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/' + task + '_Polar_AngMag.png'
plt.savefig(outfp, bbox_inches='tight')

