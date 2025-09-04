import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import os

# Load the dataset
file_path = "data/absorption_times.csv"

# Check if the file exists
if not os.path.exists(file_path):
    raise FileNotFoundError("Error: 'data/absorption_times.csv' not found.")

# Read the CSV file
df = pd.read_csv(file_path)

# Ensure required columns exist
if not {"N", "r", "AbsorptionTime"}.issubset(df.columns):
    raise ValueError("Error: Missing required columns in absorption_times.csv")

# Select a specific N and r for analysis
selected_N = 500   # Change this value as needed
selected_r = 1.2   # Change this value as needed

# Filter data for the selected parameters
subset_df = df[(df["N"] == selected_N) & (df["r"] == selected_r)]

# Check if we have enough data points
if len(subset_df) < 30:
    raise ValueError("Error: Not enough data points for analysis.")

# Extract the fixation times
fixation_times = subset_df["AbsorptionTime"].dropna()

# Compute histogram
plt.figure(figsize=(8,6))
plt.hist(fixation_times, bins=90, density=True, alpha=0.6, color='b', label="Empirical Distribution")

# Fit Gumbel distribution
gumbel_params = stats.gumbel_r.fit(fixation_times)
x = np.linspace(min(fixation_times), max(fixation_times), 100)
pdf = stats.gumbel_r.pdf(x, *gumbel_params)

# Overlay fitted Gumbel distribution
plt.plot(x, pdf, 'r-', label=f"Fitted Gumbel Distribution (params={[f'{param:.1f}' for param in gumbel_params]})")

# Labels and title
plt.xlabel("Fixation Time")
plt.ylabel("Density")
plt.title(f"Histogram of Fixation Times (N={selected_N}, r={selected_r})")
plt.legend()
plt.grid(True)

# Save the plot
plt.savefig("fixation_time_histogram.png")
plt.show()

