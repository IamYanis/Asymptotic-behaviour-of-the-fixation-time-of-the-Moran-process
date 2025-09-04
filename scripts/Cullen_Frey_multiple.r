# Load necessary packages
if (!require("fitdistrplus")) install.packages("fitdistrplus", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("moments")) install.packages("moments", dependencies = TRUE)
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)

# Load required libraries
library(fitdistrplus)
library(ggplot2)
library(moments)
library(dplyr)

# Define file path
data_file <- "data/absorption_times.csv"

# Check if file exists
if (!file.exists(data_file)) {
  stop("Error: The file 'absorption_times.csv' was not found. Please check the file path.")
}

# Read the CSV file
absorption_data <- read.csv(data_file)

# Rename columns for clarity
colnames(absorption_data) <- c("N", "r", "Sim_ID", "FixationTime")

# Remove missing values
absorption_data <- absorption_data %>% filter(!is.na(FixationTime))

# Create a folder for the plots if it doesn't exist
output_folder <- "Cullen_Frey_plot"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Get unique (N, r) combinations
unique_combinations <- unique(absorption_data[, c("N", "r")])

fixation_transform <- function(x, N, r) {
  # Example transformation: normalize by (N * r)
  return((1-r)*x +r*log(0.5) - log(N/2) - log(1-r))
}

# Loop through each unique (N, r) pair and generate Cullen & Frey plots
for (i in 1:nrow(unique_combinations)) {
  chosen_N <- unique_combinations$N[i]
  chosen_r <- unique_combinations$r[i]
  
  # Extract data for this (N, r)
  subset_data <- absorption_data %>% filter(N == chosen_N, r == chosen_r) %>% pull(FixationTime)
  
  # Check if there are enough data points
  if (length(subset_data) < 10) {
    cat("Skipping (N =", chosen_N, ", r =", chosen_r, "): Not enough data points.\n")
    next  # Skip this iteration if not enough data
  }
  
  # Apply transformation using the current N and r
  transformed_data <- fixation_transform(subset_data, chosen_N, chosen_r)

  # Define file name dynamically
  filename <- paste0(output_folder, "cullen_frey_N", chosen_N, "_r", chosen_r, ".png")
  
  # Save each plot separately
  png(filename)
  descdist(transformed_data, boot = 1000)  # Generate Cullen & Frey graph for this dataset
  
  # Overlay Gumbel distribution manually
  par(new = TRUE)  # Keep the Cullen and Frey plot
  points(1.3, 5.4 - 0.9, col = "blue", pch = 8, cex = 2)
  text(1.3, 5.4 - 0.9, labels = "Gumbel", pos = 3, col = "blue")
  
  dev.off()
  
  cat("\n Cullen and Frey graph saved as:", filename, "for (N =", chosen_N, ", r =", chosen_r, ")\n")
}
