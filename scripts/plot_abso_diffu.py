import os
import pandas as pd
import matplotlib.pyplot as plt

# --- Paths ---
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path = os.path.join(project_root, "data", "out.csv")

# --- Load data ---
df = pd.read_csv(data_path)

# Keep only absorbed trajectories
df = df[df["abs_state"] != -1]

# Separate models
df_A = df[df["model"] == "A"]
df_B = df[df["model"] == "B"]

# --- Plot ---
plt.figure(figsize=(10, 6))

plt.hist(df_A["abs_time"], bins=50, alpha=0.6, label="Model A (with drift)", density=True)
#plt.hist(df_B["abs_time"], bins=50, alpha=0.6, label="Model B (no drift)", density=True)

plt.xlabel("Absorption time")
plt.ylabel("Density")
plt.title("Histogram of Absorption Times")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()