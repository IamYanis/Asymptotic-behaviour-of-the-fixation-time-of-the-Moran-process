# Load required libraries
library(fitdistrplus)
library(ggplot2)
library(readr)
library(MASS)   # For distribution fitting
library(evd)    # For Gumbel distribution fitting
library(goftest) # For goodness-of-fit test

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
selected_r <- 1.1   # Change this value as needed

# Filter data for the selected parameters
subset_data <- abs_data[abs_data$N == selected_N & abs_data$r == selected_r, ]

# Check if we have enough data points
if (nrow(subset_data) < 30) {
  stop("Error: Not enough data points for analysis.")
}

# Extract the fixation times
fixation_times <- subset_data$AbsorptionTime

# Debugging: Check for missing values
num_missing <- sum(is.na(fixation_times))
print(paste("Number of missing values:", num_missing))

# Remove NA values
fixation_times <- na.omit(fixation_times)

# Check for extreme values
print(summary(fixation_times))
hist(fixation_times, breaks = 30, main = "Histogram of Absorption Times", col = "lightblue", border = "white")

# Ensure data is not empty after removing missing values
if (length(fixation_times) < 30) {
  stop("Error: Not enough valid data points after cleaning.")
}

# Generate Cullen and Frey plot
descdist(fixation_times, boot = 1000)

# Fit a Gumbel distribution
fit_gumbel <- fgev(fixation_times)  # From 'evd' package
loc <- fit_gumbel$estimate["location"]
scale <- fit_gumbel$estimate["scale"]

# Define Gumbel density function
dgumbel <- function(x, loc, scale) {
  z <- (x - loc) / scale
  (1/scale) * exp(-(z + exp(-z)))
}

# Define Gumbel cumulative distribution function (CDF)
pgumbel <- function(x, loc, scale) {
  exp(-exp(-(x - loc) / scale))
}

# Generate histogram and compare with Gumbel
ggplot(data.frame(fixation_times), aes(x = fixation_times)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightblue", alpha = 0.6) +
  stat_function(fun = dgumbel, args = list(loc = loc, scale = scale), color = "red", linewidth = 1.2) +
  labs(title = "Empirical vs Gumbel Distribution",
       x = "Absorption Time", y = "Density") +
  theme_minimal() +
  xlim(min(fixation_times), max(fixation_times))  # Avoid missing value warnings

# Perform Kolmogorov-Smirnov (KS) test
if (length(unique(fixation_times)) > 1) {
  ks_test <- ks.test(fixation_times, function(x) pgumbel(x, loc, scale))
  print(ks_test)
} else {
  print("Warning: KS test skipped due to insufficient unique values.")
}

# Save the plot
ggsave("empirical_vs_gumbel.png")

# Print completion message
print("Comparison with Gumbel distribution complete. Plot saved as 'empirical_vs_gumbel.png'.")
