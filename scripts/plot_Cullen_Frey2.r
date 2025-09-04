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

# Print first few rows to verify structure
print(head(absorption_data))

# Rename columns for clarity (assuming they are in the order: N, r, Sim_ID, FixationTime)
colnames(absorption_data) <- c("N", "r", "Sim_ID", "FixationTime")

# Remove missing values
absorption_data <- absorption_data %>% filter(!is.na(FixationTime))

# List unique values of N and r
unique_N <- unique(absorption_data$N)
unique_r <- unique(absorption_data$r)

# Print available options
cat("Available values for N:", unique_N, "\n")
cat("Available values for r:", unique_r, "\n")

# Select specific N and r values to analyze (modify these values as needed)
chosen_N <- unique_N[1]  # Change this to the desired N
chosen_r <- unique_r[1]  # Change this to the desired r

# Filter data for the selected N and r
filtered_data <- absorption_data %>%
  filter(N == chosen_N, r == chosen_r)

# Extract fixation times for analysis
fixation_times <- filtered_data$FixationTime

# Check if there is enough data
if (length(fixation_times) < 10) {
  stop("Not enough data points for Cullen and Frey analysis. Choose a different (N, r).")
}

# Compute empirical skewness and kurtosis
empirical_skewness <- skewness(fixation_times)
empirical_kurtosis <- kurtosis(fixation_times)

# Print empirical values
cat("Empirical Skewness for (N =", chosen_N, ", r =", chosen_r, "):", empirical_skewness, "\n")
cat("Empirical Kurtosis for (N =", chosen_N, ", r =", chosen_r, "):", empirical_kurtosis, "\n")

# Cullen and Frey Graph for the chosen (N, r)
png(paste0("cullen_frey_N", chosen_N, "_r", chosen_r, ".png"))
descdist(fixation_times, boot = 1000)
dev.off()

# Display plot in RStudio
descdist(fixation_times, boot = 1000)

cat("\n✅ Cullen and Frey graph saved as 'cullen_frey_N", chosen_N, "_r", chosen_r, ".png'\n")
