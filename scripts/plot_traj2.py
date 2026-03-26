import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Inputs you can tweak ---
file_path = "data/moran_simulation_results.csv"
N_plot = 300
gamma_plot = 0.5            # <- choose the same gamma you used in the C code
sample_trajectories = 30    # how many Moran paths to overlay
time_range = (0.0, 30.0)    # x-axis range
atol = 1e-9                 # tolerance for float comparison on r

# --- Load ---
try:
    df = pd.read_csv(file_path)
except FileNotFoundError:
    print(f"Error: File '{file_path}' not found. Run the C simulation first.")
    raise SystemExit

# Normalize column names (just in case)
df.columns = df.columns.str.strip()

# Compute r(N) for the chosen N and gamma
r_plot = 1.0 - N_plot ** (-gamma_plot)

# Filter rows for this N and (approximately) this r, and time window
mask_common = (
    (df["N"] == N_plot) &
    (np.isclose(df["r"], r_plot, atol=atol)) &
    (df["Time"] >= time_range[0]) &
    (df["Time"] <= time_range[1])
)

df_traj = df[mask_common & (df["Type"].str.strip() == "Moran")].copy()
df_ode  = df[mask_common & (df["Type"].str.strip() == "ODE")].copy()

# Sort for nicer plotting
df_traj.sort_values(by=["Sim_ID", "Time"], inplace=True)
df_ode.sort_values(by="Time", inplace=True)

# Checks
if df_traj.empty:
    print(f"Warning: No Moran data for N={N_plot}, gamma={gamma_plot} (r≈{r_plot:.6g}) in {time_range}.")
    raise SystemExit
if df_ode.empty:
    print(f"Warning: No ODE data for N={N_plot}, gamma={gamma_plot} (r≈{r_plot:.6g}) in {time_range}.")
    raise SystemExit

# Choose a subset of simulations to overlay
sim_ids = df_traj["Sim_ID"].unique()[:sample_trajectories]

plt.figure(figsize=(8, 5))

# Moran step trajectories
for i, sim_id in enumerate(sim_ids):
    sim = df_traj[df_traj["Sim_ID"] == sim_id]
    if not sim.empty:
        # label only the first for legend cleanliness
        label = "Moran process" if i == 0 else None
        plt.step(sim["Time"].to_numpy(), sim["Value"].to_numpy(), where="post", alpha=0.25, label=label)

# ODE curve
plt.plot(df_ode["Time"].to_numpy(), df_ode["Value"].to_numpy(), linestyle="--", label="Logistic ODE")

# Labels and title
plt.xlim(*time_range)
plt.xlabel("Time")
plt.ylabel("Normalized population (X/N)")
plt.title(f"Moran vs ODE (N={N_plot}, γ={gamma_plot}, r≈{r_plot:.6g})")
plt.legend()
plt.tight_layout()
plt.show()
