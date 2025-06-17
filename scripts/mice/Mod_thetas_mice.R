library(circular)
# ----- Load baseline thetas -----
files <- list.files("/scratch/users/apines/data/mouse",
                    pattern = "^m.*_1_Thetas\\.csv$",
                    full.names = TRUE)
all_data_nodrug <- lapply(files, read.csv)
nodrug_theta <- array(0, c(length(files), 503))
for (i in 1:length(files)) {
  nodrug_theta[i, ] <- all_data_nodrug[[i]]$Var1
}
nodrug_theta <- data.frame(nodrug_theta)
nodrug_theta$Drug <- 0

# ----- Load baseline rhos -----
files <- list.files("/scratch/users/apines/data/mouse",
                    pattern = "^m.*_1_Rhos\\.csv$",
                    full.names = TRUE)
all_rho_nodrug <- lapply(files, read.csv)
nodrug_rho <- array(0, c(length(files), 503))
for (i in 1:length(files)) {
  nodrug_rho[i, ] <- all_rho_nodrug[[i]]$Var1
}
nodrug_rho <- data.frame(nodrug_rho)

# ----- Load drug thetas -----
files <- list.files("/scratch/users/apines/data/mouse",
                    pattern = "^m.*_[2-6]_Thetas\\.csv$",
                    full.names = TRUE)
all_data_drug <- lapply(files, read.csv)
drug_theta <- array(0, c(length(files), 503))
for (i in 1:length(files)) {
  drug_theta[i, ] <- all_data_drug[[i]]$Var1
}
drug_theta <- data.frame(drug_theta)
drug_theta$Drug <- 1

# ----- Load drug rhos -----
files <- list.files("/scratch/users/apines/data/mouse",
                    pattern = "^m.*_[2-6]_Rhos\\.csv$",
                    full.names = TRUE)
all_rho_drug <- lapply(files, read.csv)
drug_rho <- array(0, c(length(files), 503))
for (i in 1:length(files)) {
  drug_rho[i, ] <- all_rho_drug[[i]]$Var1
}
drug_rho <- data.frame(drug_rho)

# ----- Combine everything -----
all_theta <- rbind(nodrug_theta, drug_theta)
all_rho   <- rbind(nodrug_rho,   drug_rho)

# initialize outputs
delta_x <- numeric(503)
delta_y <- numeric(503)
delta_rho <- numeric(503)
delta_theta <- numeric(503)

for (pixel in 1:503) {
  # Extract angle and rho values for this pixel
  theta_vals <- all_theta[, pixel]  # radians
  rho_vals   <- all_rho[, pixel]    # already calculated in MATLAB
  group      <- all_theta$Drug
  # Split by condition
  theta_nodrug <- theta_vals[group == 0]
  rho_nodrug   <- rho_vals[group == 0]
  theta_drug   <- theta_vals[group == 1]
  rho_drug     <- rho_vals[group == 1]
  # Number of samples
  n_nodrug <- length(theta_nodrug)
  n_drug   <- length(theta_drug)
  # Compute average resultant vector for NoDrug
  x_nodrug <- sum(rho_nodrug * cos(theta_nodrug)) / n_nodrug
  y_nodrug <- sum(rho_nodrug * sin(theta_nodrug)) / n_nodrug
  # Compute average resultant vector for Drug
  x_drug <- sum(rho_drug * cos(theta_drug)) / n_drug
  y_drug <- sum(rho_drug * sin(theta_drug)) / n_drug
  # Compute difference (Drug - NoDrug)
  dx <- x_drug - x_nodrug
  dy <- y_drug - y_nodrug
  # Save difference vector in Cartesian and Polar
  delta_x[pixel] <- dx
  delta_y[pixel] <- dy
  delta_rho[pixel] <- sqrt(dx^2 + dy^2)
  delta_theta[pixel] <- atan2(dy, dx)  # radians
}

resultant_diff <- data.frame(
  Pixel = paste0("Pixel", 1:503),
  Delta_X = delta_x,
  Delta_Y = delta_y,
  Delta_Rho = delta_rho,
  Delta_Theta = delta_theta
)

# Save to file
write.csv(resultant_diff, "~/mouse_resultant_vector_differences.csv", row.names = FALSE)
