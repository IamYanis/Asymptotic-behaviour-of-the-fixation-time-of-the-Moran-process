import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


data = pd.read_csv("data/moran_simulation_results.csv")

#extract values
N_values = sorted(data["N"].unique())
r_values = sorted(data["r"].unique())

#print("N values in dataset:", N_values)

n_sim = 10000  # This should match N_SIM in main.c

#Deviation plots
fig1, axes = plt.subplots(1, figsize=(12, 5))

scaling_exponents = {}

dev_data = data[data["Type"] == "Deviation"]  

for r in r_values:
    subset = dev_data[dev_data["r"] == r]

    N = subset["N"]
    mean_dev = subset["Value"]
    variance_dev = subset["std_sup_dev"] ** 2  

    std_err = np.sqrt(variance_dev / n_sim)  
    CI_95 = 1.96 * std_err  # 95% confidence interval

    log_N = np.log(N)
    log_mean_dev = np.log(mean_dev)

    SE_log_dev = std_err / mean_dev
    CI_log_95 = 1.96 * SE_log_dev

    slope, intercept, _, _, _ = linregress(N, log_mean_dev)
    scaling_exponents[r] = slope

    axes.errorbar(log_N, log_mean_dev, yerr=CI_log_95, marker='o', linestyle='-', capsize=5, alpha = 0.6, label=f"r={r}")

    line_styles = ['--', '-', '-.']
    style = line_styles[r_values.index(r) % len(line_styles)]

    if len(log_N) >= 4:  # we need at least 4 points to use iloc[2] and iloc[-2]
        # Use more central points to compute slope
        x_mid1, x_mid2 = log_N.iloc[2], log_N.iloc[-2]
        y_mid1, y_mid2 = log_mean_dev.iloc[2], log_mean_dev.iloc[-2]
        slope_empirical = (y_mid2 - y_mid1) / (x_mid2 - x_mid1)

        # Use slope to draw a line across the full range
        x_start, x_end = log_N.iloc[0], log_N.iloc[-1]
        y_start = y_mid1 + slope_empirical * (x_start - x_mid1)
        y_end = y_mid1 + slope_empirical * (x_end - x_mid1)

        line_styles = ['--', '-', '-.']
        style = line_styles[r_values.index(r) % len(line_styles)]
        line_label = f"slope ≈ {slope_empirical:.3f}"

        axes.plot([x_start, x_end], [y_start, y_end], linestyle=style, color='black', label=line_label)


axes.set_xlabel("Log N")
axes.set_ylabel(r"$\log E[\sup_{t \in [0,T(N)]} |X(t)/N - x(t)|]$")
#axes.set_title("Log supremum deviation")
axes.legend()
#axes.grid()

plt.tight_layout()
plt.show() 

#Print Scaling Exponents
for r, slope in scaling_exponents.items():
    print(f"For r = {r}, estimated scaling exponent: {slope:.4f}")

