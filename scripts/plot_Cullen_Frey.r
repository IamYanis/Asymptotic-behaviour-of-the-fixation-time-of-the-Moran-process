# Load required libraries
library(fitdistrplus)
library(ggplot2)
library(readr)

# Read the absorption times data
file_path <- "data/absorption_times.csv"

# Check if the file exists
if (!file.exists(file_path)) {
  stop("Error: 'data/absorption_times.csv' not found.")
}

# Load the dataset
abs_data <- read_csv(file_path)

# Ensure required columns exist
if (!all(c("N", "r", "AbsorptionTime") %in% colnames(abs_data))) {
  stop("Error: Missing required columns in absorption_times.csv")
}

# Select a specific N and r for analysis
selected_N <- 100   # Change this value as needed
selected_r <- 0.85   # Change this value as needed

# Filter data for the selected parameters
subset_data <- abs_data[abs_data$N == selected_N & abs_data$r == selected_r, ]

# Check if we have enough data points
if (nrow(subset_data) < 30) {
  stop("Error: Not enough data points for Cullen and Frey analysis")
}

# Extract the fixation times
fixation_times <- subset_data$AbsorptionTime

# Generate Cullen and Frey plot
descdist(fixation_times, boot = 1000)

# Save the plot
png("cullen_frey_plot.png")
descdist(fixation_times, boot = 1000)
dev.off()

# Print completion message
print("Cullen and Frey plot generated and saved as 'cullen_frey_plot.png'.")
