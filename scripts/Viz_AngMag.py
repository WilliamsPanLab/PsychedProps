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

# import magnitude/angle time series
MagAngL_fp='/scratch/users/apines/data/mdma/' + subj + '/' + sesh + '/' + subj + '_' + sesh + '_task-rs_p2mm_masked_Bold_and_Angles_L.mat'
MagAngR_fp='/scratch/users/apines/data/mdma/' + subj + '/' + sesh + '/' + subj + '_' + sesh + '_task-rs_p2mm_masked_Bold_and_Angles_R.mat'
MagAng_L=sio.loadmat(MagAngL_fp)['FullMatrix_L']
MagAng_R=sio.loadmat(MagAngR_fp)['FullMatrix_R']
# isolate magnitude and angles
Mag_L=MagAng_L[:,:,0]
Ang_L=MagAng_L[:,:,1]
Mag_R=MagAng_R[:,:,0]
Ang_R=MagAng_R[:,:,1]


# Define the number of angle bins and magnitude bins
num_angle_bins = 18  # 0-180 degrees divided into 18 bins (10 degrees each)
num_magnitude_bins = 8  # Adjust as needed

# Create angle bins and magnitude bins
angle_bins = np.linspace(0, np.pi, num_angle_bins + 1)  # Adjust for 0-180 degrees
# Define percentiles for magnitude bins (adjust as needed)
percentiles = np.linspace(0, 100, num_magnitude_bins + 1)  # 0-100% divided into num_magnitude_bins percentiles
# Calculate the magnitude bins based on percentiles
magnitude_bins = np.percentile(Mag_L, percentiles)

# Bin the data
binned_magnitudes = np.zeros((num_angle_bins, num_magnitude_bins))
for i in range(num_angle_bins):
    angle_min = angle_bins[i]
    angle_max = angle_bins[i + 1]
    mask = (Ang_L_rad >= angle_min) & (Ang_L_rad < angle_max)
    binned_magnitudes[i, :], _ = np.histogram(Mag_L[mask], bins=magnitude_bins)

# Create a polar contour plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

# Set the contour levels and colormap (adjust as needed)
contour = ax.contourf(np.linspace(0, np.pi, num_angle_bins), magnitude_bins, binned_magnitudes.T, cmap='viridis')

# Add a colorbar
cbar = plt.colorbar(contour, ax=ax)
cbar.set_label('Magnitude')

# Customize the plot if needed (e.g., labels, title, etc.)
ax.set_thetamin(0)  # Set the starting angle
ax.set_thetamax(180)  # Set the ending angle
ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
ax.set_theta_zero_location("N")  # Set the zero angle to the north


outfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/ACL_right_null.png'
plt.savefig(outfp,bbox_inches='tight')
