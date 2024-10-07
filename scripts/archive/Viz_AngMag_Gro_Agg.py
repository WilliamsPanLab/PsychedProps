import scipy.io as sio
import nibabel as nb
import numpy as np
import nilearn.plotting as plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import sys
import os
from PIL import Image

tasks = ['rs1', 'rs2', 'wm', 'gambling']

# Create a list to store filenames and labels for each session
sessions = ['bv', 'pl', 'm1', 'm2']
session_filenames = [[] for _ in range(len(sessions))]

for task in tasks:
    for i, session in enumerate(sessions):
        filename = f'/oak/stanford/groups/leanew1/users/apines/data/{session}_{task}_DMN_2dhist_merged.csv'
        session_filenames[i].append(filename)

# Loop over sessions
for i, session in enumerate(sessions):
    # Load data for the current session
    session_data = [np.genfromtxt(filename, delimiter=',') for filename in session_filenames[i]]
    
    # Combine binned percentile magnitudes for each task
    combined_data = np.sum(session_data, axis=0)  # Adjust combining method as needed
    
    # Create a polar contour plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    
    # Define extent (note radial axis is in radians, but labeled as degrees for interp.)
    extent = [0, np.pi, 0, 100]
    
    # Create polar plot using contourf
    contour = ax.contourf(combined_data.T, extent=extent, cmap='inferno')
    
    # Add a colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Percentile Magnitude')
    
    # Customize the plot if needed (e.g., labels, title, etc.)
    ax.set_thetamin(0)  # Set the starting angle
    ax.set_thetamax(180)  # Set the ending angle
    ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
    ax.set_theta_zero_location("N")  # Set the zero angle to the north

    # Output file path
    outfp = f'/oak/stanford/groups/leanew1/users/apines/data/p50/Group_AllTasks_Polar_AngMag_{session}.png'
    
    # Save the plot
    plt.savefig(outfp, bbox_inches='tight')
    plt.close()  # Close the figure to avoid overlapping plots

# Load data for the placebo (pl) session
data_pl = np.sum([np.genfromtxt(filename, delimiter=',') for filename in session_filenames[1]], axis=0)
total_measurements_pl = np.sum([np.genfromtxt(filename, delimiter=',') for filename in session_filenames[1]])

# Load data for the 80mg (m1) session
data_m1 = np.sum([np.genfromtxt(filename, delimiter=',') for filename in session_filenames[2]], axis=0)
total_measurements_m1 = np.sum([np.genfromtxt(filename, delimiter=',') for filename in session_filenames[2]])

# Load data for the m2 session
data_m2 = np.sum([np.genfromtxt(filename, delimiter=',') for filename in session_filenames[3]], axis=0)
total_measurements_m2 = np.sum([np.genfromtxt(filename, delimiter=',') for filename in session_filenames[3]])

# Normalize the distributions within each session
normalized_data_pl = data_pl / total_measurements_pl
normalized_data_m1 = data_m1 / total_measurements_m1
normalized_data_m2 = data_m2 / total_measurements_m2

# Calculate the differences between placebo (pl) and 80mg (m1) and between placebo (pl) and m2
difference_pl_m1 = normalized_data_pl - normalized_data_m1
difference_pl_m2 = normalized_data_pl - normalized_data_m2

# convert to percent
difference_pl_m1=difference_pl_m1*100
difference_pl_m2=difference_pl_m2*100

# Output file paths for the difference histograms
outfp_difference_pl_m1 = '/oak/stanford/groups/leanew1/users/apines/data/p50/Difference_Histogram_pl_vs_m1.png'
outfp_difference_pl_m2 = '/oak/stanford/groups/leanew1/users/apines/data/p50/Difference_Histogram_pl_vs_m2.png'

# colormap contour levels
contour_levels = np.linspace(-0.03, 0.03, 11)

# Create a 2D histogram of the differences (m1 vs. pl)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
contour_pl_m1 = ax.contourf(difference_pl_m1.T, extent=extent,levels=contour_levels, cmap='seismic')

cbar_pl_m1 = plt.colorbar(contour_pl_m1, ax=ax)
cbar_pl_m1.set_label('Difference Magnitude (pl - m2)')
ax.set_thetamin(0)  # Set the starting angle
ax.set_thetamax(180)  # Set the ending angle
ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
ax.set_theta_zero_location("N")  # Set the zero angle to the north

# Save the 2D difference histogram (m1 vs. pl)
plt.savefig(outfp_difference_pl_m1, bbox_inches='tight')
plt.close()  # Close the figure to avoid overlapping plots

# Create a 2D histogram of the differences (m2 vs. pl)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
contour_pl_m2 = ax.contourf(difference_pl_m2.T, extent=extent,levels=contour_levels, cmap='seismic')
cbar_pl_m2 = plt.colorbar(contour_pl_m2, ax=ax)
cbar_pl_m2.set_label('Difference Magnitude (pl - m2)')
ax.set_thetamin(0)  # Set the starting angle
ax.set_thetamax(180)  # Set the ending angle
ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
ax.set_theta_zero_location("N")  # Set the zero angle to the north

# Save the 2D difference histogram (m2 vs. pl)
plt.savefig(outfp_difference_pl_m2, bbox_inches='tight')
plt.close()  # Close the figure to avoid overlapping plots


