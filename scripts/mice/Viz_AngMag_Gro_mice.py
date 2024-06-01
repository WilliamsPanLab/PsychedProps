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


# load in group-level angle/magnitude 2-d histograms
fn1='/oak/stanford/groups/leanew1/users/apines/data/mice_1_DMN_2dhist_merged.csv'
fn2='/oak/stanford/groups/leanew1/users/apines/data/mice_2_DMN_2dhist_merged.csv'
fn3='/oak/stanford/groups/leanew1/users/apines/data/mice_3_DMN_2dhist_merged.csv'
fn4='/oak/stanford/groups/leanew1/users/apines/data/mice_4_DMN_2dhist_merged.csv'
fn5='/oak/stanford/groups/leanew1/users/apines/data/mice_5_DMN_2dhist_merged.csv'
fn6='/oak/stanford/groups/leanew1/users/apines/data/mice_6_DMN_2dhist_merged.csv'

# Define the number of angle bins and magnitude bins
num_angle_bins = 10  # 0-180 degrees divided into 18 bins (10 degrees each)
num_magnitude_bins = 10

# Create angle bins and magnitude bins
angle_bins = np.linspace(0, np.pi, num_angle_bins+1)  # Adjust for 0-180 degrees
# Define percentiles for magnitude bins (adjust as needed)
percentiles = np.linspace(0, 100, num_magnitude_bins+1)  # 0-100% divided into num_magnitude_bins percentiles

# run 1
binned_percentile_magnitudes=np.genfromtxt(fn1,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/mice_LSD_Polar_AngMag_1.png'
plt.savefig(outfp, bbox_inches='tight')


# run 2
binned_percentile_magnitudes=np.genfromtxt(fn2,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/mice_LSD_Polar_AngMag_2.png'
plt.savefig(outfp, bbox_inches='tight')


# run 3
binned_percentile_magnitudes=np.genfromtxt(fn3,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/mice_LSD_Polar_AngMag_3.png'
plt.savefig(outfp, bbox_inches='tight')


# run 4
binned_percentile_magnitudes=np.genfromtxt(fn4,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/mice_LSD_Polar_AngMag_4.png'
plt.savefig(outfp, bbox_inches='tight')


# run 5
binned_percentile_magnitudes=np.genfromtxt(fn5,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/mice_LSD_Polar_AngMag_5.png'
plt.savefig(outfp, bbox_inches='tight')

# run 6
binned_percentile_magnitudes=np.genfromtxt(fn6,delimiter=',')
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
outfp = '/oak/stanford/groups/leanew1/users/apines/data/p50/mice_LSD_Polar_AngMag_6.png'
plt.savefig(outfp, bbox_inches='tight')
