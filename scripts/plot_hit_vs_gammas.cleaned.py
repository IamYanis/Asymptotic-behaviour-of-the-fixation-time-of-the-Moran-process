import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --------- choose the population size N you want to inspect ---------
target_N = 1000

# --------- load (robust path + robust parsing) ---------
here = Path(__file__).resolve().parent
data_path = (here / ".." / "data" / "absorption_times.cleaned.csv").resolve()  # <-- garantit le bon fichier
df = pd.read_csv(data_path, low_memory=False)

# --------- normalize column names ---------
df.columns = [c.strip() for c in df.columns]
if "N" not in df.columns:
    raise KeyError(f"Colonne 'N' absente. Colonnes présentes: {df.columns.tolist()}")

# --------- coerce types (critique pour le filtrage) ---------
df["N_norm"] = pd.to_numeric(df["N"], errors="coerce").round().astype("Int64")
df["AbsorptionTime"] = pd.to_numeric(df.get("AbsorptionTime", np.nan), errors="coerce")
df["gamma"] = pd.to_numeric(df.get("gamma", np.nan), errors="coerce")
df["AbsorbingState"] = pd.to_numeric(df.get("AbsorbingState", np.nan), errors="coerce")

# --------- quick debug: affiche les N disponibles ---------
avail_N = np.sort(df["N_norm"].dropna().unique())
print("Available N values:", avail_N[:25], "... total:", len(avail_N))

# --------- filter ---------
dfN = df[df["N_norm"] == target_N].copy()
if dfN.empty:
    raise ValueError(
        f"No rows found for N={target_N}. N disponibles (premiers 25): {avail_N[:25]}.\n"
        f"Vérifie que tu lis bien: {data_path}"
    )

# (optionnel) garde les temps finis
dfN = dfN[np.isfinite(dfN["AbsorptionTime"])]

# --------- aggregate counts by gamma ---------
cnt0  = dfN[dfN["AbsorbingState"] == 0].groupby("gamma").size()
cntN  = dfN[dfN["AbsorbingState"] == target_N].groupby("gamma").size()
cntNA = dfN[dfN["AbsorbingState"] == -1].groupby("gamma").size()

counts = pd.DataFrame({
    "Hits at 0": cnt0,
    f"Hits at N={target_N}": cntN,
    "Not absorbed": cntNA
}).fillna(0).astype(int).sort_index()

# --------- plot ---------
gammas = counts.index.to_numpy(dtype=float)
cnt_0  = counts["Hits at 0"].to_numpy()
cnt_N  = counts[f"Hits at N={target_N}"].to_numpy()
cnt_NA = counts["Not absorbed"].to_numpy()

fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(gammas)); w = 0.25
ax.bar(x - w, cnt_0, width=w, label="Hits at 0 (extinction)")
ax.bar(x,     cnt_N, width=w, label=f"Hits at N={target_N} (fixation)")
ax.bar(x + w, cnt_NA, width=w, label="Not absorbed")
ax.set_xlabel("γ"); ax.set_ylabel("Number of runs")
ax.set_title(f"Outcome counts vs γ for N={target_N}")
ax.set_xticks(x); ax.set_xticklabels([f"{g:.3f}" for g in gammas], rotation=45, ha="right")
ax.legend()
plt.tight_layout()
plt.show()
