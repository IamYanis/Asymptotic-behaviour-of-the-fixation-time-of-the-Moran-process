import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
N_selected = 50000
r_selected = 0.8
csv_path = "data/absorption_times.csv"

# Force correct column interpretation (ignore broken header)
cols = ["N", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState", "Absorbed"]
df = pd.read_csv(csv_path, header=None, names=cols, skiprows=1)

# Convert to numeric
for c in cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Keep valid absorption times
df = df[np.isfinite(df["AbsorptionTime"])]

# OPTIONAL: keep only absorbed runs
df = df[df["Absorbed"] == 1]

# Select N and r (float-safe)
mask = (df["N"] == N_selected) & np.isclose(df["r"], r_selected)
times = df.loc[mask, "AbsorptionTime"].to_numpy()

if len(times) == 0:
    print("Existing N values:", sorted(df["N"].unique()))
    print("Existing r values:", sorted(df["r"].unique()))
    raise ValueError("No matching data found.")

# Plot histogram (density)
plt.hist(times, bins=100, density=True, edgecolor="black", alpha=0.7)
plt.xlabel("Absorption Time")
plt.ylabel("Density")
plt.title(f"N={N_selected}, r={r_selected}")
plt.tight_layout()
plt.show()
