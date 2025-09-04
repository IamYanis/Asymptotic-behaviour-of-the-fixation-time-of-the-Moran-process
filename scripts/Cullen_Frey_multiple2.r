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

# Select a specific value of r (so we only compare across N)
chosen_r <- 1.1  # Change this to the desired r value

# Filter data for the selected r
filtered_data <- absorption_data %>% filter(r == chosen_r)

# Compute skewness and kurtosis for each N (only for selected r)
summary_data <- filtered_data %>%
  group_by(N) %>%
  summarise(
    Skewness = skewness(FixationTime),
    Kurtosis = kurtosis(FixationTime),
    .groups = 'drop'
  )

# Compute skewness squared (needed for Cullen and Frey graph)
summary_data$Skewness2 <- summary_data$Skewness^2

# Define theoretical distribution locations
theoretical_distributions <- data.frame(
  Name = c("Normal", "Uniform", "Exponential", "Gamma", "Lognormal", "Beta", "Gumbel"),
  Skewness2 = c(0, 0, 2, 3, 4, 0.5, 1.3),  # Square of skewness values
  Kurtosis = c(3, 1.8, 9, 8, 7, 2.5, 5.4)
)

# Plot Cullen and Frey graph for all (N) values on one graph
ggplot() +
  # Add theoretical distribution points
  geom_point(data = theoretical_distributions, aes(x = Skewness2, y = Kurtosis), 
             color = "black", size = 3, shape = 8) +  # Star marker for theoretical points
  geom_text(data = theoretical_distributions, aes(x = Skewness2, y = Kurtosis, label = Name), 
            vjust = -0.5, hjust = 1.1, size = 4, color = "black") +
  
  # Add empirical data points for all (N, r) on the same graph
  geom_point(data = summary_data, aes(x = Skewness2, y = Kurtosis, color = factor(N)), size = 4) +
  geom_text(data = summary_data, aes(x = Skewness2, y = Kurtosis, label = paste("N=", N)), 
            vjust = -0.5, hjust = 0.5, size = 3.5) +

  # Labels and styling
  labs(title = paste("Cullen and Frey Graph for r =", chosen_r),
       x = "Square of Skewness",
       y = "Kurtosis",
       color = "N (Population Size)") +
  theme_minimal()

# Save the plot
ggsave(paste0("Cullen_Frey_plot/cullen_frey_overlay_r", chosen_r, ".png"), width = 8, height = 6)

cat("\n✅ Cullen and Frey graph with multiple N values saved as 'Cullen_Frey_plot/cullen_frey_overlay_r", chosen_r, ".png'\n")
