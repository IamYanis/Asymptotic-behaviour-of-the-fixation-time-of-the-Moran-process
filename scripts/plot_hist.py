import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----------------- Load Absorption Times Data ----------------- #
abs_data = pd.read_csv("data/absorption_times.csv")

# Debug: Check for invalid values
print("Checking for invalid values in AbsorptionTime column...")
print(abs_data["AbsorptionTime"].describe())
print(abs_data["AbsorptionTime"].isna().sum(), "NaN values found")
print(np.isinf(abs_data["AbsorptionTime"]).sum(), "Inf values found")

# Remove NaN and Inf values before plotting
abs_data = abs_data[np.isfinite(abs_data["AbsorptionTime"])]

# Select 4 values of N for the histograms
selected_N = [50, 100, 300, 10000]

# ----------------- Plot 2x2 Histograms ----------------- #
fig, hist_axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)

r_values = sorted(abs_data["r"].unique())  # Extract unique r values

for idx, N in enumerate(selected_N):
    row, col = divmod(idx, 2)  # Convert index to 2D subplot grid
    ax = hist_axes[row, col]  # Select correct subplot

    subset = abs_data[abs_data["N"] == N]
    
    #for r in r_values[:5]:  # Plot histograms for first 4 r-values
    for r in r_values:
        times = subset[subset["r"] == r]["AbsorptionTime"]

        # Remove invalid values
        times = times[np.isfinite(times)]

        ax.hist(times, bins=120, alpha=0.5, histtype='stepfilled', label=f"r={r:.2f}", edgecolor="black", range=(0, 500))

    ax.set_xlabel("Absorption Time")
    ax.set_title(f"N = {N}")
    ax.legend()

# Set shared labels and formatting
hist_axes[0, 0].set_ylabel("Frequency")
hist_axes[1, 0].set_ylabel("Frequency")

#plt.suptitle("Histograms of Absorption Times for Different N")
plt.tight_layout()
plt.show()
