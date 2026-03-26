import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gumbel_r

FILE_PATH = "data/absorption_times.csv"

def load_absorption_csv(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Error: '{path}' not found.")

    # Peek at the file to detect header/column count reliably
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        raise ValueError("Error: absorption_times.csv is empty.")

    first = lines[0]
    has_header = first.lower().startswith("n,")

    # Find the most common number of comma-separated fields in the data lines
    data_lines = lines[1:] if has_header else lines
    if not data_lines:
        raise ValueError("Error: absorption_times.csv has a header but no data rows.")

    counts = [len(ln.split(",")) for ln in data_lines]
    col_count = max(set(counts), key=counts.count)

    if col_count == 6:
        # Matches current main.c fprintf: N,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed
        names = ["N", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState", "Absorbed"]
        df = pd.read_csv(path, header=None, names=names, skiprows=1 if has_header else 0)
    elif col_count == 7:
        # Matches header: N,gamma,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed
        df = pd.read_csv(path)
    else:
        raise ValueError(
            f"Unsupported column count ({col_count}). Expected 6 or 7 columns. "
            f"Saw counts like: {sorted(set(counts))}"
        )

    # Coerce numeric columns safely
    for c in ["N", "r", "AbsorptionTime"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["N", "r", "AbsorptionTime"])

    return df

def make_grid_axes(k: int, base=4.0, sharey=True):
    """
    Create a near-square grid of subplots for k panels.
    Returns (fig, axes_flat, nrows, ncols).
    """
    ncols = math.ceil(math.sqrt(k))
    nrows = math.ceil(k / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(base * ncols, base * nrows), sharey=sharey)
    axes = np.array(axes).ravel()
    return fig, axes, nrows, ncols

def hide_unused_axes(axes, k, nrows, ncols):
    for j in range(k, nrows * ncols):
        axes[j].set_visible(False)

def label_left_column(axes, ncols, ylabel: str):
    for ax in axes[::ncols]:
        if ax.get_visible():
            ax.set_ylabel(ylabel)

# ============================================================
# LOAD DATA
# ============================================================
abs_data = load_absorption_csv(FILE_PATH)

# Optional: keep only absorbed trajectories if present
if "Absorbed" in abs_data.columns:
    abs_data["Absorbed"] = pd.to_numeric(abs_data["Absorbed"], errors="coerce")
    abs_data = abs_data[abs_data["Absorbed"] == 1]

# Pick N values that actually exist in the data
selected_N = [100, 400, 1000, 5000]
available_N = sorted(abs_data["N"].astype(int).unique().tolist())
selected_N = [N for N in selected_N if N in available_N]
if not selected_N:
    raise ValueError(f"No selected_N values found in file. Available N: {available_N}")

r_values = sorted(abs_data["r"].unique())

L_range = 5
prec = 10

# ============================================================
# SECOND FIGURE: RAW ABSORPTION-TIME DENSITY HISTOGRAMS (SQUARE GRID)
# ============================================================

def g_bd_extinction(t, r, X0, N):
    a = (1.0 / r - 1.0)  # >0 when r<1
    e = np.exp(-a * t)
    num = (1.0 - e) ** (N - X0 - 1)
    den = (r - e) ** (N - X0 + 1)
    return (N - X0) * a ** 2 * e * (num / den)

def f_hit_m(t, lbda, mu, N):
    a = (lbda - mu)
    b = mu / lbda
    e = np.exp(-a * t)
    num = (1 - e) ** (N - 2)
    den = (1 - b * e) ** (N + 1)
    return (1 - b) ** 2 * num / den * N * e

def Q1_asymptotic(t, lbda, mu, N):
    if lbda <= mu:
        raise ValueError("This approximation requires lbda > mu (supercritical).")

    t = np.asarray(t, dtype=float)

    alpha = lbda - mu
    beta = mu / lbda

    # centering
    tN = (np.log(N) + np.log(1.0 - beta)) / alpha

    # scaled variable
    y = alpha * (t - tN)

    return (1.0 - beta) * alpha * np.exp(-y) * np.exp(-np.exp(-y))

def f_gumbel_shift_scale(t, r, X0_bis, N):
    c = -r * np.log(1.0 - X0_bis / N) + np.log(X0_bis * (1.0 - r))
    return (1.0 - r) * np.exp(-((1.0 - r) * t - c + np.exp((r - 1.0) * t + c)))

def fft_convolution_density(r, X0, X0_bis, N, tmax=50.0, dt=1e-2):
    t = np.arange(0.0, tmax - 60 + dt, dt)

    g = Q1_asymptotic(t, 1.0 / r, 1, N=N - X0_bis)
    f = f_gumbel_shift_scale(t, r=r, X0_bis=X0_bis, N=N)

    L = len(t)
    nfft = 1 << (2 * L - 1).bit_length()

    G = np.fft.rfft(g, n=nfft)
    F = np.fft.rfft(f, n=nfft)
    h = np.fft.irfft(G * F, n=nfft)[:L] * dt

    h = np.maximum(h, 0.0)

    area = np.trapz(h, t)
    if area > 0:
        h /= area

    return t, h, g, f


k2 = len(selected_N)
fig2, axes2, nrows2, ncols2 = make_grid_axes(k2, base=4.0, sharey=True)

for i, N in enumerate(selected_N):
    ax = axes2[i]
    subsetN = abs_data[abs_data["N"].astype(int) == int(N)]

    #X0 = N - 1
    X0 = N/2

    for r in r_values:
        times = subsetN.loc[np.isclose(subsetN["r"], r), "AbsorptionTime"].to_numpy()
        if times.size == 0:
            continue


        ax.hist(
            times,
            bins=50,
            density=True,
            alpha=0.5,
            histtype="stepfilled",
            edgecolor="black",
            label=f"r={r:.2f}",
        )

    ax.set_title(f"N = {N}")
    ax.set_xlabel("Absorption time")
    ax.legend(fontsize=8)

label_left_column(axes2, ncols2, "Density")
hide_unused_axes(axes2, k2, nrows2, ncols2)
plt.tight_layout()
plt.show()

# ============================================================
# DATA QUALITY REPORT
# ============================================================

def quality_report(df: pd.DataFrame, time_max_rule="N"):
    """
    time_max_rule:
      - "N"  -> flags AbsorptionTime > N (per row)
      - number -> flags AbsorptionTime > that number
    """
    print("\n================ DATA QUALITY REPORT ================\n")
    print("Rows:", len(df))
    print("Columns:", list(df.columns))

    cols_to_check = [c for c in ["N", "r", "AbsorptionTime", "AbsorbingState", "Absorbed"] if c in df.columns]
    print("\nNaN counts:")
    print(df[cols_to_check].isna().sum())

    if "AbsorbingState" in df.columns:
        st = pd.to_numeric(df["AbsorbingState"], errors="coerce")
        df2 = df.copy()
        df2["AbsorbingState_num"] = st

        t0 = (df2["AbsorbingState_num"] == 0).sum()
        tn = (df2["AbsorbingState_num"] == df2["N"]).sum() if "N" in df2.columns else 0
        other = (~df2["AbsorbingState_num"].isin([0]) & (df2["AbsorbingState_num"] != df2["N"])).sum()

        print("\nAbsorption destination counts (from AbsorbingState):")
        print(f"  T0 (hit 0): {t0}")
        print(f"  TN (hit N): {tn}")
        print(f"  TN+T0= {t0 + tn}")
        print(f"  Other / weird states: {other}")

    df_num = df.copy()
    df_num["N"] = pd.to_numeric(df_num["N"], errors="coerce")
    df_num["AbsorptionTime"] = pd.to_numeric(df_num["AbsorptionTime"], errors="coerce")

    if time_max_rule == "N":
        bad = df_num[
            np.isfinite(df_num["AbsorptionTime"]) & np.isfinite(df_num["N"]) &
            (df_num["AbsorptionTime"] > df_num["N"])
        ]
        rule_str = "AbsorptionTime > N"
    else:
        bad = df_num[
            np.isfinite(df_num["AbsorptionTime"]) &
            (df_num["AbsorptionTime"] > float(time_max_rule))
        ]
        rule_str = f"AbsorptionTime > {time_max_rule}"

    print(f"\nSuspiciously large times ({rule_str}): {len(bad)} rows")
    if len(bad) > 0:
        print("Top 10 largest times:")
        show_cols = [c for c in ["N", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState"] if c in bad.columns]
        print(bad.sort_values("AbsorptionTime", ascending=False)[show_cols].head(10).to_string(index=False))

    print("\nBy (N, r) summary:")
    grp = df_num.groupby(["N", "r"], dropna=False)
    summary = grp["AbsorptionTime"].agg(["count", "min", "median", "max"])
    print(summary.to_string())

    print("\n=====================================================\n")

quality_report(abs_data, time_max_rule="N")