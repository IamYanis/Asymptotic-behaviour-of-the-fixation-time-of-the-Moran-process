import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = "data/moran_simulation_results.csv"
try:
    df = pd.read_csv(file_path)
except FileNotFoundError:
    print(f"Error: File '{file_path}' not found. Run the C simulation first.")
    exit()

# Select specific (N, r) for visualization
N_plot = 20
r_plot = 0.8
time_range = (0, 200)  # Set the desired time range

# Normalize column names (remove unwanted spaces)
df.columns = df.columns.str.strip()

# Filter data for Moran process and ODE solution within the time range
df_traj = df[
    (df["N"] == N_plot) & 
    (df["r"] == r_plot) & 
    (df["Type"].str.strip() == "Moran") & 
    (df["Time"] >= time_range[0]) & 
    (df["Time"] <= time_range[1])
].copy()

df_ode = df[
    (df["N"] == N_plot) & 
    (df["r"] == r_plot) & 
    (df["Type"].str.strip() == "ODE") & 
    (df["Time"] >= time_range[0]) & 
    (df["Time"] <= time_range[1])
].copy()

# Check if data is available
if df_traj.empty:
    print(f"Warning: No Moran process data found for N={N_plot}, r={r_plot} in the time range {time_range}. Check your C simulation output.")

if df_ode.empty:
    print(f"Warning: No ODE solution found for N={N_plot}, r={r_plot} in the time range {time_range}. Check your C simulation output.")

if df_traj.empty or df_ode.empty:
    exit()

# **Assign Simulation ID based on resets in Time**
df_traj["Sim_ID"] = (df_traj["Time"].diff() < 0).cumsum()

# Get first few unique simulations
simulations = df_traj["Sim_ID"].unique()[:200]

plt.figure(figsize=(8, 5))

# Determine final time in range
final_time = time_range[1]

# Track which trajectories should be red
for i, sim_id in enumerate(simulations):
    sim_data = df_traj[df_traj["Sim_ID"] == sim_id].sort_values(by="Time")
    
    if not sim_data.empty:
        last_time = sim_data["Time"].max()
        last_value = sim_data["Value"].iloc[-1]
        
        # Plot full trajectory in black
        label = "Moran Process" if i == 0 else None
        plt.step(sim_data["Time"], sim_data["Value"], alpha=0.1, color="black", label=label)
        
        # Highlight last point in red
        plt.scatter(last_time, last_value, color="red", s=10)

# Plot ODE solution (continuous curve)
plt.plot(df_ode["Time"], df_ode["Value"], color="red", linestyle="dashed", label="Logistic ODE")

# Labels and Title
plt.xlabel("Time")
plt.ylabel("Normalized Population (X/N)")
plt.title(f"Moran Process vs Logistic ODE (N={N_plot}, r={r_plot}, Time {time_range[0]}-{time_range[1]})")
plt.legend()
plt.show()
