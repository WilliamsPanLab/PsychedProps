# Define the directory where your data files are located
data_dir <- "/scratch/users/apines/data/mdma"

# Initialize a 4x17 matrix with zeros
matrix <- matrix(0, nrow = 4, ncol = 17)

# List all subject directories
subject_dirs <- list.files(data_dir, full.names = TRUE)

# Iterate over subject directories
for (subject_index in 1:length(subject_dirs)) {
  # Extract the subject name from the directory
  subject_name <- basename(subject_dirs[subject_index])
  
  # List all session directories for the subject
  session_dirs <- list.files(subject_dirs[subject_index], full.names = TRUE)
  
  # Iterate over session directories
  for (session_index in 1:length(session_dirs)) {
    # Extract the session name from the directory
    session_name <- basename(session_dirs[session_index])
    
    # Construct the path to the text file
    txt_file <- file.path(session_dirs[session_index], paste(subject_name, session_name, "task-rs_ValidSegments_Trunc.txt", sep = "_"))
    
    # Check if the file exists
    if (file.exists(txt_file)) {
      # Read the data from the text file
      data <- read.table(txt_file, header = FALSE, sep = ",")
      
      # Add the fourth column to the corresponding position in the matrix
      matrix[session_index, subject_index] <- sum(data$V2)
    }
  }
}

# Print the resulting matrix
print(matrix)

