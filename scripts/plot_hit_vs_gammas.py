import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --------- choose the population size N you want to inspect ---------
target_N = 5000   # <-- change this to any N present in your CSV

# --------- load & filter ---------
df = pd.read_csv("data/absorption_times.csv")

# keep only this N
dfN = df[df["N"] == target_N].copy()
if dfN.empty:
    raise ValueError(f"No rows found for N={target_N}.")

# guard against non-finite times (optional)
if "AbsorptionTime" in dfN.columns:
    dfN = dfN[np.isfinite(dfN["AbsorptionTime"])]

# --------- aggregate counts by gamma ---------
cnt0  = dfN[dfN["AbsorbingState"] == 0].groupby("gamma").size()
cntN  = dfN[dfN["AbsorbingState"] == target_N].groupby("gamma").size()
cntNA = dfN[dfN["AbsorbingState"] == -1].groupby("gamma").size()

# join everything into one DataFrame
counts = pd.DataFrame({
    "Hits at 0": cnt0,
    f"Hits at N={target_N}": cntN,
    "Not absorbed": cntNA
}).fillna(0).astype(int).sort_index()

# values for plotting
gammas = counts.index.to_numpy(dtype=float)
cnt_0  = counts["Hits at 0"].to_numpy()
cnt_N  = counts[f"Hits at N={target_N}"].to_numpy()
cnt_NA = counts["Not absorbed"].to_numpy()

# --------- plot: grouped bars ---------
fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(gammas))
w = 0.25

for gamma in gammas:
    print(gamma)


ax.bar(x - w, cnt_0, width=w, label="Hits at 0 (extinction)")
ax.bar(x,     cnt_N, width=w, label=f"Hits at N={target_N} (fixation)")
ax.bar(x + w, cnt_NA, width=w, label="Not absorbed")

# cosmetics
ax.set_xlabel("γ")
ax.set_ylabel("Number of runs")
ax.set_title(f"hitting times vs γ for N={target_N}")
ax.set_xticks(x)
ax.set_xticklabels([f"{g:.3f}" for g in gammas])
ax.legend()

plt.tight_layout()
plt.show()
