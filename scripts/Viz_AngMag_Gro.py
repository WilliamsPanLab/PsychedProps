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

task = sys.argv[1]

# load in group-level angle/magnitude 2-d histograms


bvfn = '/oak/stanford/groups/leanew1/users/apines/data/bv_' + task + '_DMN_2dhist_merged.csv'
plfn = '/oak/stanford/groups/leanew1/users/apines/data/pl_' + task + '_DMN_2dhist_merged.csv'
m1fn = '/oak/stanford/groups/leanew1/users/apines/data/m1_' + task + '_DMN_2dhist_merged.csv'
m2fn = '/oak/stanford/groups/leanew1/users/apines/data/m2_' + task + '_DMN_2dhist_merged.csv'

# Define the number of angle bins and magnitude bins
num_angle_bins = 10  # 0-180 degrees divided into 18 bins (10 degrees each)
num_magnitude_bins = 10

# Create angle bins and magnitude bins
angle_bins = np.linspace(0, np.pi, num_angle_bins+1)  # Adjust for 0-180 degrees
# Define percentiles for magnitude bins (adjust as needed)
percentiles = np.linspace(0, 100, num_magnitude_bins+1)  # 0-100% divided into num_magnitude_bins percentiles

# bv
binned_percentile_magnitudes=np.genfromtxt(bvfn,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/Group_' + task + '_Polar_AngMag_bv.png'
plt.savefig(outfp, bbox_inches='tight')

# pl
binned_percentile_magnitudes=np.genfromtxt(plfn,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/Group_' + task + '_Polar_AngMag_pl.png'
plt.savefig(outfp, bbox_inches='tight')

# m1
binned_percentile_magnitudes=np.genfromtxt(m1fn,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/Group_' + task + '_Polar_AngMag_m1.png'
plt.savefig(outfp, bbox_inches='tight')

# m2
binned_percentile_magnitudes=np.genfromtxt(m2fn,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/Group_' + task + '_Polar_AngMag_m2.png'
plt.savefig(outfp, bbox_inches='tight')

