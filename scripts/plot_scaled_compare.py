import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ABS_PATH = "data/absorption_times.csv"
OUT_PATH = "data/out.csv"          # limit diffusion absorption times
FIG_DIR = "figures"

# -------------------------------
# Load absorption_times.csv (your robust version)
# -------------------------------
def load_absorption_csv(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Error: '{path}' not found.")

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if not lines:
        raise ValueError("Error: absorption_times.csv is empty.")

    first = lines[0]
    has_header = first.lower().startswith("n,")

    data_lines = lines[1:] if has_header else lines
    if not data_lines:
        raise ValueError("Error: absorption_times.csv has a header but no data rows.")

    counts = [len(ln.split(",")) for ln in data_lines]
    col_count = max(set(counts), key=counts.count)

    if col_count == 6:
        names = ["N", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState", "Absorbed"]
        df = pd.read_csv(path, header=None, names=names, skiprows=1 if has_header else 0)
    elif col_count == 7:
        df = pd.read_csv(path)  # N,gamma,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed
    else:
        raise ValueError(f"Unsupported column count ({col_count}). Expected 6 or 7.")

    for c in ["N", "r", "AbsorptionTime"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["N", "r", "AbsorptionTime"]).copy()

    if "Absorbed" in df.columns:
        df["Absorbed"] = pd.to_numeric(df["Absorbed"], errors="coerce")
        df = df[df["Absorbed"] == 1].copy()

    df["N"] = df["N"].astype(int)
    return df


# -------------------------------
# Load out.csv (limit diffusion)
# -------------------------------
def load_limit_csv(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Error: '{path}' not found. Put your limit diffusion file at {OUT_PATH}."
        )

    df = pd.read_csv(path)

    # Try to find the time column
    # Preferred: AbsorptionTime, else time, else t, else last numeric col
    if "AbsorptionTime" in df.columns:
        tcol = "AbsorptionTime"
    elif "time" in df.columns:
        tcol = "time"
    elif "t" in df.columns:
        tcol = "t"
    else:
        tmp = df.copy()
        for c in tmp.columns:
            tmp[c] = pd.to_numeric(tmp[c], errors="ignore")
        num_cols = [c for c in tmp.columns if pd.api.types.is_numeric_dtype(tmp[c])]
        if not num_cols:
            raise ValueError(f"Could not infer a numeric time column in {path}. Columns={list(df.columns)}")
        tcol = num_cols[-1]

    df = df.copy()
    df[tcol] = pd.to_numeric(df[tcol], errors="coerce")
    df = df.dropna(subset=[tcol]).rename(columns={tcol: "AbsorptionTime"}).copy()

    # Optional columns (if present)
    for c in ["N", "r"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Optional absorbed flag
    if "Absorbed" in df.columns:
        df["Absorbed"] = pd.to_numeric(df["Absorbed"], errors="coerce")
        df = df[df["Absorbed"] == 1].copy()

    if "N" in df.columns:
        df["N"] = df["N"].astype("Int64")

    return df


# -------------------------------
# Plot helpers
# -------------------------------
def make_grid_axes(k: int, base=4.0, sharey=True):
    ncols = math.ceil(math.sqrt(k))
    nrows = math.ceil(k / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(base * ncols, base * nrows), sharey=sharey)
    return fig, np.array(axes).ravel(), nrows, ncols

def hide_unused_axes(axes, k, nrows, ncols):
    for j in range(k, nrows * ncols):
        axes[j].set_visible(False)

def label_left_column(axes, ncols, ylabel: str):
    for ax in axes[::ncols]:
        if ax.get_visible():
            ax.set_ylabel(ylabel)


# ============================================================
# MAIN
# ============================================================
def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    abs_df = load_absorption_csv(ABS_PATH)
    lim_df = load_limit_csv(OUT_PATH)

    selected_N = [1000, 5000, 10000, 50000]
    available_N = sorted(abs_df["N"].unique().tolist())
    selected_N = [N for N in selected_N if N in available_N]
    if not selected_N:
        raise ValueError(f"No selected N found. Available N: {available_N}")

    r_values = sorted(abs_df["r"].unique())

    fig, axes, nrows, ncols = make_grid_axes(len(selected_N), base=4.2, sharey=True)



    for i, N in enumerate(selected_N):
        ax = axes[i]
        sub_abs_N = abs_df[abs_df["N"] == int(N)]

        # Filter limit data by N if it has an N column; otherwise use all limit times
        sub_lim_N = lim_df
        if "N" in lim_df.columns and lim_df["N"].notna().any():
            sub_lim_N = lim_df[lim_df["N"].astype("Int64") == int(N)]

        for r in r_values:
            # Empirical times for this (N,r)
            emp = sub_abs_N.loc[np.isclose(sub_abs_N["r"], r), "AbsorptionTime"].to_numpy()
            if emp.size == 0:
                continue

            # Limit times: filter by r if present; otherwise treat as same limit for all r
            sub_lim_Nr = sub_lim_N
            if "r" in sub_lim_N.columns:
                sub_lim_Nr = sub_lim_N.loc[np.isclose(sub_lim_N["r"], r)]

            lim = sub_lim_Nr["AbsorptionTime"].to_numpy()

            # Use the SAME bins for both (union support)
            if lim.size > 0:
                combined = np.concatenate([emp, lim])
            else:
                combined = emp

            # robust bins from combined data
            bins = np.histogram_bin_edges(combined, bins=50)

            # Empirical histogram (filled)
            ax.hist(emp, bins=bins, 
                    density=True, 
                    alpha=0.45, histtype="stepfilled", edgecolor="black", label=f"Empirical r={r:.2f}")

            # Limit histogram (outline)
            if lim.size > 0:
                ax.hist(lim, bins=bins,
                         density=True, 
                         histtype="step", linewidth=2, label=f"Limit r={r:.2f}" if "r" in sub_lim_N.columns else "Limit" )
                #None
            else:
                # If you want, comment this out to keep plots cleaner
                pass

        ax.set_title(f"N = {N}")
        ax.set_xlabel("Absorption time (scaled)")
        ax.legend(fontsize=8)

    label_left_column(axes, ncols, "Density")
    hide_unused_axes(axes, len(selected_N), nrows, ncols)
    plt.tight_layout()

    out_fig = os.path.join(FIG_DIR, "hist_empirical_vs_limit.png")
    plt.savefig(out_fig, dpi=200)
    print(f"Saved: {out_fig}")
    plt.show()


if __name__ == "__main__":
    main()