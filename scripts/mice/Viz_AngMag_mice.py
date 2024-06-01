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
MagAng_fp='/scratch/users/apines/p50_mice/proc/20200228/' + subj + '_' + sesh + '_masked_Bold_and_Angles.mat'
MagAng=sio.loadmat(MagAng_fp)['FullMatrix']
# split out rows where no angular optical flow value was recorded (most posterior row of the brain)
rows_all_180 = np.all(MagAng[:,:,1] == 180, axis=1)
rows_valid = rows_all_180==False
# isolate magnitude and angles
Mag=MagAng[rows_valid,:,0]
Ang=MagAng[rows_valid,:,1]
Ang_rad = np.radians(Ang)

# Define the number of angle bins and magnitude bins
num_angle_bins = 10  # 0-180 degrees divided into 18 bins (10 degrees each)
num_magnitude_bins = 10

# Create angle bins and magnitude bins
angle_bins = np.linspace(0, np.pi, num_angle_bins+1)  # Adjust for 0-180 degrees
# Define percentiles for magnitude bins (adjust as needed)
percentiles = np.linspace(0, 100, num_magnitude_bins+1)  # 0-100% divided into num_magnitude_bins percentiles
magnitude_percentiles = np.percentile(Mag, percentiles)

# Bin the data using global percentiles for each angle bin: right
binned_percentile_magnitudes = np.zeros((num_angle_bins, num_magnitude_bins))
for i in range(num_angle_bins):
    angle_min = angle_bins[i]
    angle_max = angle_bins[i + 1]
    mask = (Ang_rad >= angle_min) & (Ang_rad < angle_max)
    # Calculate percentile bins for the magnitude data using global percentiles
    bin_edges = magnitude_percentiles
    # Create histograms using the calculated bin edges
    binned_percentile_magnitudes[i, :], _ = np.histogram(Mag[mask], bins=bin_edges)


# save out for group-level merging
output_csv_path='/scratch/users/apines/p50_mice/proc/20200228/' + subj + '_' + sesh +  '_DMN_2dhist.csv'
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
outfp='/scratch/users/apines/p50_mice/proc/20200228/' + subj + '_' + sesh + '_Polar_AngMag.png'
plt.savefig(outfp, bbox_inches='tight')

