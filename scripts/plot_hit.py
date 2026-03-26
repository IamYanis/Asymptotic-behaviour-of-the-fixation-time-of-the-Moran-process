# This code shows the density of the two hitting times T_0 and T_N




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

# ---------- Load & clean ----------
df = pd.read_csv("data/absorption_times.csv")

# Keep only absorbed runs (we're plotting hitting times)
if "Absorbed" in df.columns:
    df = df[df["Absorbed"] == 1]

# Drop non-finite times
df = df[np.isfinite(df["AbsorptionTime"])]

if df.empty:
    raise ValueError("No absorbed runs with finite AbsorptionTime found.")

# ---------- Group by (N, r, gamma) ----------
groups = list(df.groupby(["N", "r", "gamma"], sort=True))
n_groups = len(groups)

# Layout
ncols = 2 if n_groups > 1 else 1
nrows = math.ceil(n_groups / ncols)

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 4 * nrows), squeeze=False)
axes = axes.flatten()

for ax, ((N, r, gamma), g) in zip(axes, groups):
    # Split by absorbing state
    times0 = g.loc[g["AbsorbingState"] == 0, "AbsorptionTime"].to_numpy()
    timesN = g.loc[g["AbsorbingState"] == N, "AbsorptionTime"].to_numpy()

    # Skip empty panel
    if len(times0) == 0 and len(timesN) == 0:
        ax.set_visible(False)
        continue

    # Shared bin edges
    all_times = np.concatenate([times for times in (times0, timesN) if len(times) > 0])
    lo = 0.0
    hi = float(np.percentile(all_times, 99.5))
    if not np.isfinite(hi) or hi <= lo:
        hi = float(all_times.max())

    bins = min(80, max(10, int(np.sqrt(len(all_times)))))

    if len(times0) > 0:
        ax.hist(times0, bins=bins, range=(lo, hi), density=True, alpha=0.5,
                label=f"Extinction at 0 (n={len(times0)})")
    if len(timesN) > 0:
        ax.hist(timesN, bins=bins, range=(lo, hi), density=True, alpha=0.5,
                label=f"Fixation at N (n={len(timesN)})")

    ax.set_title(f"N={N}, r={r:.6f}, γ={gamma:.3f}")
    #ax.set_xlabel("Hitting time")
    ax.set_ylabel("Density")
    ax.legend()

# Hide unused axes
for ax in axes[n_groups:]:
    ax.set_visible(False)


plt.tight_layout()
plt.show()
