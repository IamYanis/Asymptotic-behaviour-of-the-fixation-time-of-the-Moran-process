import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

# ----------------------------
# PARAMETERS
# ----------------------------
N_selected = 1000
i0 = 1                      # initial mutant count in theory (change if your sims start elsewhere)
csv_path = "data/absorption_times.csv"
bins = 300
gamma_slice = slice(4, None)  # [4:] ; use slice(2,5) for [2:5]
xmax = 1500

# ----------------------------
# THEORY: Continuous-time Moran absorption-time density
# ----------------------------
def moran_Q(N: int, r: float) -> np.ndarray:
    """
    Subgenerator Q on transient states {1,...,N-1} for continuous-time Moran process.
    Fitness: mutant=r, wildtype=1.
    Time scale c=1.
    """
    m = N - 1
    Q = np.zeros((m, m), dtype=float)

    for i in range(1, N):  # i = 1..N-1
        denom = r * i + (N - i)
        lam = (r * i / denom) * ((N - i) / N)   # i -> i+1
        mu  = ((N - i) / denom) * (i / N)       # i -> i-1

        idx = i - 1
        Q[idx, idx] = -(lam + mu)
        if i < N - 1:
            Q[idx, idx + 1] = lam
        if i > 1:
            Q[idx, idx - 1] = mu

    return Q

def absorption_pdf(N: int, r: float, i0: int, tgrid: np.ndarray) -> np.ndarray:
    """
    Exact PDF of absorption time T for CTMC, starting at i0 mutants,
    absorbing at {0, N}.
    f(t) = alpha^T exp(Qt) (-Q 1)
    """
    Q = moran_Q(N, r)
    m = N - 1
    one = np.ones(m)

    alpha = np.zeros(m)
    alpha[i0 - 1] = 1.0

    exit_rates = -Q @ one  # (-Q 1)

    pdf = np.empty_like(tgrid, dtype=float)
    for k, t in enumerate(tgrid):
        v = alpha @ expm(Q * t)     # (m,)
        pdf[k] = float(v @ exit_rates)
    return pdf

# ----------------------------
# DATA: Read your simulation CSV
# ----------------------------
cols = ["N", "gamma", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState", "Absorbed"]
df = pd.read_csv(csv_path, header=None, names=cols, skiprows=1)

for c in cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df[np.isfinite(df["AbsorptionTime"])]
df = df[df["Absorbed"] == 1]

dfN = df[df["N"] == N_selected].copy()
if dfN.empty:
    print("Existing N values:", sorted(df["N"].dropna().unique()))
    raise ValueError(f"No data found for N={N_selected}")

gamma_values = np.sort(dfN["gamma"].dropna().unique())[gamma_slice]
print("gamma values found:", gamma_values)

# Common bins for fair comparison (for histograms)
t_all = dfN["AbsorptionTime"].to_numpy()
t_all = t_all[np.isfinite(t_all)]
bin_edges = np.linspace(t_all.min(), t_all.max(), bins + 1)

# Colormap for vivid colors
cmap = plt.cm.plasma
colors = cmap(np.linspace(0, 1, len(gamma_values)))

# ----------------------------
# PLOT: histogram + theoretical curve per gamma
# ----------------------------
plt.figure(figsize=(9, 6))

for gamma, color in zip(gamma_values, colors):
    mask = np.isclose(dfN["gamma"], gamma)
    times = dfN.loc[mask, "AbsorptionTime"].to_numpy()
    if times.size < 5:
        continue

    # r corresponding to this gamma (should be constant within gamma)
    r_vals = dfN.loc[mask, "r"].dropna().unique()
    if r_vals.size == 0:
        continue
    r = float(r_vals[0])

    # Histogram (simulation)
    plt.hist(
        times,
        bins=bin_edges,
        density=True,
        alpha=0.35,
        color=tuple(color),
        label=f"sim γ={gamma:g} (r={r:.4g})"
    )

    # Theory curve (exact density)
    # Build a t-grid on the visible range for overlay
    tmax = min(xmax, np.max(times))
    tgrid = np.linspace(0, tmax, 400)
    pdf = absorption_pdf(N_selected, r, i0, tgrid)

    plt.plot(
        tgrid,
        pdf,
        linewidth=2.0,
        color=tuple(color),
        linestyle="--",
        label=f"theory γ={gamma:g}"
    )

plt.xlim(0, xmax)
plt.xlabel("Absorption Time")
plt.ylabel("Density")
plt.title(f"Absorption-time density: simulation vs exact theory (N={N_selected}, i0={i0})")
plt.legend(ncol=2, fontsize=9)
plt.tight_layout()
plt.show()
