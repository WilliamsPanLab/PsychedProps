import scipy.io as sio
import nibabel as nb
import numpy as np
import nilearn.plotting as plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import sys
import os
from PIL import Image


# load in aggregated data
data_1=np.genfromtxt('/oak/stanford/groups/leanew1/users/apines/data/mice_1_DMN_2dhist_merged.csv',delimiter=',')
data_2=np.genfromtxt('/oak/stanford/groups/leanew1/users/apines/data/mice_2_DMN_2dhist_merged.csv',delimiter=',')
data_3=np.genfromtxt('/oak/stanford/groups/leanew1/users/apines/data/mice_3_DMN_2dhist_merged.csv',delimiter=',')
data_4=np.genfromtxt('/oak/stanford/groups/leanew1/users/apines/data/mice_4_DMN_2dhist_merged.csv',delimiter=',')
data_5=np.genfromtxt('/oak/stanford/groups/leanew1/users/apines/data/mice_5_DMN_2dhist_merged.csv',delimiter=',')
data_6=np.genfromtxt('/oak/stanford/groups/leanew1/users/apines/data/mice_6_DMN_2dhist_merged.csv',delimiter=',')

# sep into pre and post drug
predrug=data_1
total_measurements_predrug=np.sum(predrug)

# post drug
postdrug=data_2+data_3+data_4+data_5+data_6
total_measurements_postdrug=np.sum(postdrug)

# Normalize the distributions within each session
normalized_data_pre = predrug / total_measurements_predrug
normalized_data_post = postdrug / total_measurements_postdrug

# Calculate the differences between placebo (pl) and 80mg (m1) and between placebo (pl) and m2
difference_prepost = normalized_data_pre - normalized_data_post

# convert to percent
difference_prepost=difference_prepost*100

# Output file paths for the difference histograms
outfp_difference_prepost = '/oak/stanford/groups/leanew1/users/apines/data/p50/Difference_Histogram_sobermice_minus_acidmice.png'

# colormap contour levels
contour_levels = np.linspace(-0.24, 0.24, 11)
extent = [0, np.pi, 0, 100]

# Create a 2D histogram of the differences (m1 vs. pl)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
contour_pl_m1 = ax.contourf(difference_prepost.T, extent=extent,levels=contour_levels, cmap='seismic')

cbar_pl_m1 = plt.colorbar(contour_pl_m1, ax=ax)
cbar_pl_m1.set_label('Difference Magnitude (pre - post)')
ax.set_thetamin(0)  # Set the starting angle
ax.set_thetamax(180)  # Set the ending angle
ax.set_theta_direction(-1)  # Rotate clockwise (adjust as needed)
ax.set_theta_zero_location("N")  # Set the zero angle to the north

# Save the 2D difference histogram (m1 vs. pl)
plt.savefig(outfp_difference_prepost, bbox_inches='tight')
plt.close()  # Close the figure to avoid overlapping plots



