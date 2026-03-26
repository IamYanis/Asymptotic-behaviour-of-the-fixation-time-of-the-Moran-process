"Plot the histogram of the density of the absorption time for different value of r = 1+\delta/N^gamma for N fixed"


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
N_selected = 1000
csv_path = "data/absorption_times.csv"
bins = 300

# Read (ignore broken header)
cols = ["N", "gamma", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState", "Absorbed"]
df = pd.read_csv(csv_path, header=None, names=cols, skiprows=1)

# Convert to numeric
for c in cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Keep valid times and absorbed runs
df = df[np.isfinite(df["AbsorptionTime"])]
df = df[df["Absorbed"] == 1]

# Filter by N
dfN = df[df["N"] == N_selected].copy()

if dfN.empty:
    print("Existing N values:", sorted(df["N"].dropna().unique()))
    raise ValueError(f"No data found for N={N_selected}")

# r values available for this N
gamma_values = np.sort(dfN["gamma"].dropna().unique())
#gamma_values = np.sort(dfN["gamma"].dropna().unique())[4:]
print("r (gamma) values found:", gamma_values)

# Use common bins for fair comparison
t_all = dfN["AbsorptionTime"].to_numpy()
bin_edges = np.linspace(t_all.min(), t_all.max(), bins + 1)

plt.figure(figsize=(9, 6))

cmap = plt.cm.plasma   # very saturated
colors = cmap(np.linspace(0, 1, len(gamma_values)))


for gamma, color in zip(gamma_values, colors):
    times = dfN.loc[np.isclose(dfN["gamma"], gamma), "AbsorptionTime"].to_numpy()
    if times.size < 5:
        continue
    plt.hist(times, bins=bin_edges, density=True, histtype="stepfilled", alpha=0.5,  color = color, edgecolor = "black", label=f"gamma={gamma:g}")

plt.xlim(0, 1500)
plt.xlabel("Absorption Time")
plt.ylabel("Density")
plt.title(f"Absorption time density by r (N={N_selected})")
plt.legend()
plt.tight_layout()
plt.show()
