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

# Select a specific value of r
chosen_r <- 1.1  # Change this to your desired r value

# Filter data for the selected r
filtered_data <- absorption_data %>% filter(r == chosen_r)

# Get unique N values
unique_N <- unique(filtered_data$N)

# Set up plot
png("cullen_frey_multiple_N.png")

# Generate the Cullen and Frey graph for the first N (as a base)
first_N <- unique_N[1]
first_data <- filtered_data %>% filter(N == first_N) %>% pull(FixationTime)
descdist(first_data, boot = 1000)

# Loop through remaining N values and overlay points
par(new = TRUE)  # Keep the base plot

for (N_value in unique_N) {
  # Extract data for the current N
  fixation_times <- filtered_data %>% filter(N == N_value) %>% pull(FixationTime)

  # Compute skewness and kurtosis manually
  skew_val <- skewness(fixation_times)
  kurt_val <- kurtosis(fixation_times)
  skew_sq <- skew_val^2  # Square of skewness

  # Add empirical point on Cullen and Frey graph
  points(skew_sq, kurt_val, col = "red", pch = 19)  # Red points for each N
  text(skew_sq, kurt_val, labels = paste("N=", N_value), pos = 3, cex = 0.8)
}

dev.off()

cat("\n✅ Cullen and Frey graph for multiple N values (fixed r =", chosen_r, ") saved as 'cullen_frey_multiple_N.png'\n")
