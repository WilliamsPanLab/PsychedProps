import scipy.io as sio
import nibabel as nb
import numpy as np
import nilearn.plotting as plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib
matplotlib.use('Agg')
import sys
import os
from PIL import Image

# Create a list to store filenames and labels for each session
sessions = ['pre_psil', 'during_psil', 'during_meth']
session_filenames = [[] for _ in range(len(sessions))]

for i, session in enumerate(sessions):
    filename = f'/oak/stanford/groups/leanew1/users/apines/data/{session}_DMN_2dhist_merged.mat'
    session_filenames[i].append(filename)

# Loop over sessions
for i, session in enumerate(sessions):
    # Load data for the current session
    fn=session_filenames[i]
    fnstr=fn[0]
    session_data = sio.loadmat(fnstr.strip("[]"))
    matrix_name = list(session_data.keys())[-1] 
    session_data = session_data[matrix_name]
    # Combine binned percentile magnitudes for each task
    combined_data = np.sum(session_data, axis=2)  # Adjust combining method as needed
    total_measurements = np.sum(session_data)
    normalized_data = combined_data / total_measurements
    # Create a polar contour plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    # Define extent (note radial axis is in radians, but labeled as degrees for interp.)
    extent = [0, np.pi, 0, 100]
    # Create polar plot using contourf
    contour = ax.contourf(normalized_data.T, extent=extent, cmap='inferno',vmin= 0.00945,vmax=0.01050)
    # Add a colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Percentile Magnitude')
    # Customize the plot if needed (e.g., labels, title, etc.)
    ax.set_thetamin(0)  # Set the starting angle
    ax.set_thetamax(180)  # Set the ending angle
    ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
    ax.set_theta_zero_location("N")  # Set the zero angle to the north
    # Output file path
    outfp = f'/oak/stanford/groups/leanew1/users/apines/data/p50/Group_rs_Polar_AngMag_{matrix_name}.png'
    # Save the plot
    plt.savefig(outfp, bbox_inches='tight')
    plt.close()  # Close the figure to avoid overlapping plots

# Load data for the pre (pr) session
pre_mat=sio.loadmat('/oak/stanford/groups/leanew1/users/apines/data/pre_psil_DMN_2dhist_merged.mat')
data_pr=np.sum(pre_mat['outDF_pre'], axis=2)
total_measurements_pr = np.sum(pre_mat['outDF_pre'])

# Load data for psil session
psil_mat=sio.loadmat('/oak/stanford/groups/leanew1/users/apines/data/during_psil_DMN_2dhist_merged.mat')
data_psil= np.sum(psil_mat['outDF_psil'], axis=2)
total_measurements_psil = np.sum(psil_mat['outDF_psil'])

# Load data for the methyl
meth_mat=sio.loadmat('/oak/stanford/groups/leanew1/users/apines/data/during_meth_DMN_2dhist_merged.mat')
data_meth = np.sum(meth_mat['outDF_meth'], axis=2)
total_measurements_meth = np.sum(meth_mat['outDF_meth'])

# Normalize the distributions within each session
normalized_data_pr = data_pr / total_measurements_pr
normalized_data_psil = data_psil / total_measurements_psil
normalized_data_meth = data_meth / total_measurements_meth

# Calculate the differences between placebo (pl) and 80mg (m1) and between placebo (pl) and m2
difference_pr_psil = normalized_data_pr - normalized_data_psil
difference_pr_meth = normalized_data_pr - normalized_data_meth
# convert to percent
difference_pr_psil=difference_pr_psil*100
difference_pr_meth=difference_pr_meth*100

# Output file paths for the difference histograms
outfp_difference_pr_psil = '/oak/stanford/groups/leanew1/users/apines/data/p50/Difference_Histogram_pr_vs_psil.png'
outfp_difference_pr_meth = '/oak/stanford/groups/leanew1/users/apines/data/p50/Difference_Histogram_pr_vs_meth.png'

# adding merged version for clarity
outfp_difference_merged = '/oak/stanford/groups/leanew1/users/apines/data/p50/Difference_Histogram_pl_vs_drug.png'

# colormap contour levels
contour_levels = np.linspace(-0.025, 0.025, 11)

# Create a 2D histogram of the differences (m1 vs. pl)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
contour_pl_m1 = ax.contourf(difference_pr_psil.T, extent=extent,levels=contour_levels, cmap='seismic')

cbar_pl_m1 = plt.colorbar(contour_pl_m1, ax=ax)
cbar_pl_m1.set_label('Difference Magnitude (pre - psil)')
ax.set_thetamin(0)  # Set the starting angle
ax.set_thetamax(180)  # Set the ending angle
ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
ax.set_theta_zero_location("N")  # Set the zero angle to the north

# Save the 2D difference histogram (m1 vs. pl)
plt.savefig(outfp_difference_pr_psil, bbox_inches='tight')
plt.close()  # Close the figure to avoid overlapping plots

# Create a 2D histogram of the differences (m2 vs. pl)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
contour_pl_m2 = ax.contourf(difference_pr_meth.T, extent=extent,levels=contour_levels, cmap='seismic')
cbar_pl_m2 = plt.colorbar(contour_pl_m2, ax=ax)
cbar_pl_m2.set_label('Difference Magnitude (pre - meth)')
ax.set_thetamin(0)  # Set the starting angle
ax.set_thetamax(180)  # Set the ending angle
ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
ax.set_theta_zero_location("N")  # Set the zero angle to the north

# Save the 2D difference histogram (m2 vs. pl)
plt.savefig(outfp_difference_pr_meth, bbox_inches='tight')
plt.close()  # Close the figure to avoid overlapping plots

