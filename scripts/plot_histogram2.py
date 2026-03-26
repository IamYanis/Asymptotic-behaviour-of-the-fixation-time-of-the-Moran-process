import os
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
    # mode
    col_count = max(set(counts), key=counts.count)

    if col_count == 6:
        # Matches current main.c fprintf: N,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed
        names = ["N", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState", "Absorbed"]
        df = pd.read_csv(path, header=None, names=names, skiprows=1 if has_header else 0)
    elif col_count == 7:
        # Matches the header you intended: N,gamma,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed
        # If file truly has 7 cols, standard read works.
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported column count ({col_count}). "
                         f"Expected 6 or 7 columns. Saw counts like: {sorted(set(counts))}")

    # Coerce numeric columns safely
    for c in ["N", "r", "AbsorptionTime"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["N", "r", "AbsorptionTime"])

    return df

abs_data = load_absorption_csv(FILE_PATH)

# Optional: keep only absorbed trajectories if present
if "Absorbed" in abs_data.columns:
    abs_data["Absorbed"] = pd.to_numeric(abs_data["Absorbed"], errors="coerce")
    abs_data = abs_data[abs_data["Absorbed"] == 1]

# Pick N values that actually exist in the data (your main.c uses 10,30,100,1000)
selected_N = [400, 1000, 10000]
available_N = sorted(abs_data["N"].astype(int).unique().tolist())
selected_N = [N for N in selected_N if N in available_N]
if not selected_N:
    raise ValueError(f"No selected_N values found in file. Available N: {available_N}")

r_values = sorted(abs_data["r"].unique())

L_range = 5
prec = 10




fig, axes = plt.subplots(1, len(selected_N), figsize=(4 * len(selected_N), 4), sharey=True)
if len(selected_N) == 1:
    axes = [axes]

for i, N in enumerate(selected_N):
    subsetN = abs_data[abs_data["N"].astype(int) == int(N)]

    #X0 = np.floor(.7*N)
    #X0 = N/2
    #X0 = np.floor(N - np.sqrt(N))
    #X0 = np.floor(N - N**.25)
    X0 = N-15

    for r in r_values:
        times = subsetN.loc[np.isclose(subsetN["r"], r), "AbsorptionTime"].to_numpy()
        if times.size == 0:
            continue

        # Your original scaling (valid for r < 1, and you said r=0.8)
        transformed = (1 - r) * times + r * np.log(1 - X0/N) - np.log(X0) - np.log(1 - r)

        axes[i].hist(
            transformed,
            bins=L_range * prec,
            density=True,
            alpha=0.5,
            histtype="stepfilled",
            edgecolor="black",
            range=(-L_range, L_range + 2.5),
            label=f"r={r:.2f}",
        )

    axes[i].set_title(f"N = {N}")
    axes[i].legend(fontsize=8)

    # Standard Gumbel overlay
    x = np.linspace(-L_range, L_range, prec * 10)
    y = gumbel_r.pdf(x, loc=0, scale=1)
    axes[i].plot(x, y, linestyle="dashed", color="red", alpha=0.8, label="Std Gumbel")

axes[0].set_ylabel("Density")
plt.tight_layout()
#plt.show()

# ============================================================
# SECOND FIGURE: RAW ABSORPTION-TIME DENSITY HISTOGRAMS
# ============================================================

fig2, axes2 = plt.subplots(
    1, len(selected_N),
    figsize=(4 * len(selected_N), 4),
    sharey=True
)

if len(selected_N) == 1:
    axes2 = [axes2]

for i, N in enumerate(selected_N):
    subsetN = abs_data[abs_data["N"].astype(int) == int(N)]

    #X0 = np.floor(.7*N)
    #X0 = N/2
    #X0 = np.floor(N - np.sqrt(N))
    #X0 = np.floor(N - N**.25)
    X0 = N-15


    for r in r_values:
        times = subsetN.loc[np.isclose(subsetN["r"], r), "AbsorptionTime"].to_numpy()
        if times.size == 0:
            continue

        L_range2 = 100
        x2 = np.linspace(0, L_range2, prec * 10)
        c = - r*np.log(1 - X0/N) + np.log(X0 * (1-r))
        #c = r*np.log(N/(N-X0)) + np.log(X0 * (1-r))
        y2 = (1-r) * np.exp( - ((1-r)*x2 - c + np.exp((r-1)*x2+c)))
        #X0 = N/2
        z2 = (N-X0) * (1/r-1)**2 * np.exp(-(1/r-1)*x2) * ( (1-np.exp(-(1/r-1)*x2))**(N-X0-1) / (1/r - np.exp(-(1/r-1)*x2))**(N-X0+1))
        #pN = (r**(N-X0) - r**N)/(1-r**N)
        #R = 1/r
        #pN = (R**(X0) - R**N)/(1-R**N)
        pN_approx = r**(N-X0)
        p0 = (1-r**(N-X0))/(1-r**N)
        #p0 = (1-R**(X0))/(1-R**N)
        p0_approx = 1 - r**(N-X0)
        f2 = z2 + p0*y2
        print("r=", r, "p0 =",p0)
        #f2 = (r**(N-X0) - r**N)*z2 + (1-r**(N-X0))*y2

        #axes2[i].plot(x2, y2, linestyle="dashed", color="blue", alpha=0.8, label="T density")
        #axes2[i].plot(x2, z2, linestyle="dashed", color="blue", alpha=0.8, label="T density")
        axes2[i].plot(x2, f2, linestyle="dashed", color="blue", alpha=0.8, label="T density")
        #axes2[i].plot(x2, y2, linestyle="dashed", color="r", alpha=0.8, label="T density")


        axes2[i].hist(
            times,
            bins=50,                 # adjust if needed
            density=True,
            alpha=0.5,
            histtype="stepfilled",
            edgecolor="black",
            label=f"r={r:.2f}",
        )

    axes2[i].set_title(f"N = {N}")
    axes2[i].set_xlabel("Absorption time")
    axes2[i].legend(fontsize=8)

    

axes2[0].set_ylabel("Density")
plt.tight_layout()
plt.show()

def quality_report(df: pd.DataFrame, time_max_rule="N"):
    """
    time_max_rule:
      - "N"  -> flags AbsorptionTime > N (per row)
      - number -> flags AbsorptionTime > that number
    """
    print("\n================ DATA QUALITY REPORT ================\n")
    print("Rows:", len(df))
    print("Columns:", list(df.columns))

    # Missing values counts
    cols_to_check = [c for c in ["N", "r", "AbsorptionTime", "AbsorbingState", "Absorbed"] if c in df.columns]
    print("\nNaN counts:")
    print(df[cols_to_check].isna().sum())

    # Count T0 vs TN if AbsorbingState exists
    if "AbsorbingState" in df.columns:
        # In your C code, this should be 0 or N (sometimes stored as int)
        st = pd.to_numeric(df["AbsorbingState"], errors="coerce")
        df = df.copy()
        df["AbsorbingState_num"] = st

        t0 = (df["AbsorbingState_num"] == 0).sum()
        tn = (df["AbsorbingState_num"] == df["N"]).sum() if "N" in df.columns else 0
        other = (~df["AbsorbingState_num"].isin([0]) & (df["AbsorbingState_num"] != df["N"])).sum()

        print("\nAbsorption destination counts (from AbsorbingState):")
        print(f"  T0 (hit 0): {t0}")
        print(f"  TN (hit N): {tn}")
        print(f" TN+T0= {t0+tn}")
        print(f"  Other / weird states: {other}")

    # Flag huge times
    df_num = df.copy()
    df_num["N"] = pd.to_numeric(df_num["N"], errors="coerce")
    df_num["AbsorptionTime"] = pd.to_numeric(df_num["AbsorptionTime"], errors="coerce")

    if time_max_rule == "N":
        bad = df_num[np.isfinite(df_num["AbsorptionTime"]) & np.isfinite(df_num["N"]) &
                     (df_num["AbsorptionTime"] > df_num["N"])]
        rule_str = "AbsorptionTime > N"
    else:
        bad = df_num[np.isfinite(df_num["AbsorptionTime"]) & (df_num["AbsorptionTime"] > float(time_max_rule))]
        rule_str = f"AbsorptionTime > {time_max_rule}"

    print(f"\nSuspiciously large times ({rule_str}): {len(bad)} rows")
    if len(bad) > 0:
        print("Top 10 largest times:")
        show_cols = [c for c in ["N", "r", "Sim_ID", "AbsorptionTime", "AbsorbingState"] if c in bad.columns]
        print(bad.sort_values("AbsorptionTime", ascending=False)[show_cols].head(10).to_string(index=False))

    # Summary by (N, r): counts + max time + NaNs
    print("\nBy (N, r) summary:")
    grp = df_num.groupby(["N", "r"], dropna=False)
    summary = grp["AbsorptionTime"].agg(["count", "min", "median", "max"])
    print(summary.to_string())

    print("\n=====================================================\n")


# Call it:
quality_report(abs_data, time_max_rule="N")   # or time_max_rule=1e6 etc.
