import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----------------- PARAMETERS ----------------- #
N_selected = 10000    # choose N
r_selected = 0.8     # choose r
csv_path = "data/absorption_times.csv"

# ----------------- LOAD DATA ----------------- #
# main.c writes:
# N,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed
cols = ["N", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState", "Absorbed"]

df = pd.read_csv(csv_path, names=cols, header=0)

# Ensure numeric types
df["N"] = pd.to_numeric(df["N"], errors="coerce")
df["r"] = pd.to_numeric(df["r"], errors="coerce")
df["AbsorptionTime"] = pd.to_numeric(df["AbsorptionTime"], errors="coerce")

# Drop invalid values
df = df[np.isfinite(df["AbsorptionTime"])]

# Keep only absorbed trajectories (recommended)
df = df[df["Absorbed"] == 1]

# ----------------- SELECT DATA ----------------- #
times = df.loc[
    (df["N"] == N_selected) & (df["r"] == r_selected),
    "AbsorptionTime"
].to_numpy()

if len(times) == 0:
    raise ValueError(f"No data found for N={N_selected}, r={r_selected}")

# ----------------- PLOT HISTOGRAM ----------------- #
plt.figure(figsize=(6, 4))
plt.hist(
    times,
    bins=100,
    density=True,
    histtype="stepfilled",
    edgecolor="black",
    alpha=0.7
)

plt.xlabel("Absorption Time")
plt.ylabel("Density")
plt.title(f"Absorption Time Histogram (N={N_selected}, r={r_selected})")

plt.tight_layout()
plt.show()
