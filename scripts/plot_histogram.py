import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from scipy.stats import gumbel_r

file_path = "data/absorption_times.csv"
if not os.path.exists(file_path):
    print(f"Error: '{file_path}' not found.")
    exit()

abs_data = pd.read_csv(file_path)

if "N" not in abs_data.columns or "r" not in abs_data.columns or "AbsorptionTime" not in abs_data.columns:
    print("Error: Missing required columns in absorption_times.csv")
    exit()

#selected_N = [20, 70, 300, 1000]
selected_N = [10, 30, 100, 1000]
r_values = sorted(abs_data["r"].unique())

fig, axes = plt.subplots(1, len(selected_N), figsize=(15, 5), sharey=True)

L_range = 5
prec = 10

for i, N in enumerate(selected_N):
    subset = abs_data[abs_data["N"] == N]
    for r in r_values:
        times = subset[subset["r"] == r]["AbsorptionTime"]
        if times.empty:
            continue
        axes[i].hist( (1-r)*times +r*np.log(0.5) - np.log(N/2) - np.log(1-r), bins=L_range*prec, density=True, alpha=0.5, histtype='stepfilled', label=f"r={r:.2f}", edgecolor="black", range=(-L_range, L_range+2.5))

    #axes[i].set_xlabel("Absorption Time")
    axes[i].set_title(f"N = {N}")
    axes[i].legend()

    #loc, scale = gumbel_r.fit(times)
    loc, scale = 0, 1 #Standard Gumbel
    x = np.linspace(-L_range, L_range, prec*10)
    y = gumbel_r.pdf(x, loc=loc, scale=scale)
    axes[i].plot(x, y, linestyle="dashed", color="red", alpha=0.8, label=f"Gumbel r={r:.2f}")

    

axes[0].set_ylabel("Frequency")
#plt.suptitle("Histograms of Absorption Times")
plt.tight_layout()
plt.show()
