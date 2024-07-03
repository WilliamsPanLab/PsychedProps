import os

# Define the directory path
directory = "/scratch/users/apines/p50_mice/proc/20200310"

# Create a new list to store modified lines
new_lines = []

# Iterate through each folder in the specified directory
for folder_name in os.listdir(directory):
    # Checking if it's a folder and contains "postLSD"
    if os.path.isdir(os.path.join(directory, folder_name)) and "Ket25mgkg" in folder_name:
        # Since we are in the directory "20200228", we can directly use this date
        date = "20200310"
        new_folder_name = f"{date} {folder_name}"
        new_lines.append(new_folder_name)

# Write the modified lines to a new text file
with open("/scratch/users/apines/p50_mice/demo/ipynb/sess_ket25.txt", "w") as file:
    for line in new_lines:
        file.write(line + "\n")

