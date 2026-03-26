#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter


def robust_read_csv(path: str) -> pd.DataFrame:
    """
    Read a potentially messy CSV. Tries the fast parser first,
    then falls back to the python engine with autodetected sep,
    and finally skips bad lines if needed.
    """
    try:
        return pd.read_csv(path)
    except Exception:
        try:
            return pd.read_csv(path, engine="python", sep=None)
        except Exception:
            return pd.read_csv(path, engine="python", sep=None, on_bad_lines="skip")


def main(csv_path: str, target_N: int | None):
    # --------- load & basic checks ---------
    df = robust_read_csv(csv_path)

    required = {"N", "gamma", "AbsorbingState"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}. "
            f"Found columns: {sorted(df.columns.tolist())}"
        )

    # Auto-pick the most common N if none provided
    if target_N is None:
        target_N = Counter(df["N"]).most_common(1)[0][0]

    # keep only this N
    dfN = df[df["N"] == target_N].copy()
    if dfN.empty:
        avail = sorted(pd.unique(df["N"].dropna()).tolist())
        raise ValueError(
            f"No rows found for N={target_N}. "
            f"Available N values (first 30): {avail[:30]}{'...' if len(avail) > 30 else ''}"
        )

    # guard against non-finite times (optional)
    if "AbsorptionTime" in dfN.columns:
        dfN = dfN[np.isfinite(dfN["AbsorptionTime"])]

    # --------- aggregate counts by gamma ---------
    cnt0 = dfN[dfN["AbsorbingState"] == 0].groupby("gamma").size()
    cntN = dfN[dfN["AbsorbingState"] == target_N].groupby("gamma").size()
    cntNA = dfN[dfN["AbsorbingState"] == -1].groupby("gamma").size()

    # join everything into one DataFrame
    counts = (
        pd.DataFrame(
            {
                "Hits at 0": cnt0,
                f"Hits at N={target_N}": cntN,
                "Not absorbed": cntNA,
            }
        )
        .fillna(0)
        .astype(int)
        .sort_index()
    )

    # print a quick preview in the terminal
    print("\nAggregated counts (first 20 rows):")
    print(counts.head(20).to_string())

    # --------- plot: grouped bars ---------
    gammas = counts.index.to_numpy(dtype=float)
    cnt_0 = counts["Hits at 0"].to_numpy()
    cnt_N = counts[f"Hits at N={target_N}"].to_numpy()
    cnt_NA = counts["Not absorbed"].to_numpy()

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(gammas))
    w = 0.25

    ax.bar(x - w, cnt_0, width=w, label="Hits at 0 (extinction)")
    ax.bar(x, cnt_N, width=w, label=f"Hits at N={target_N} (fixation)")
    ax.bar(x + w, cnt_NA, width=w, label="Not absorbed")

    # cosmetics
    ax.set_xlabel("γ")
    ax.set_ylabel("Number of runs")
    ax.set_title(f"Outcome counts vs γ for N={target_N}")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g:.3f}" for g in gammas], rotation=0)
    ax.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot absorption outcome counts by γ for a chosen population size N."
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default="data/absorption_times.csv",
        help="Path to absorption_times CSV (default: data/absorption_times.csv)",
    )
    parser.add_argument(
        "--N",
        type=int,
        default=None,
        help="Population size N to filter. If omitted, uses the most common N in the CSV.",
    )
    args = parser.parse_args()
    main(args.csv, args.N)
