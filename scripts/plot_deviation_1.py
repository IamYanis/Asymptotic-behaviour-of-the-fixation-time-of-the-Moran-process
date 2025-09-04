import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


data = pd.read_csv("data/moran_simulation_results.csv")

#extract values
N_values = sorted(data["N"].unique())
r_values = sorted(data["r"].unique())

#print("N values in dataset:", N_values)

n_sim = 1000  # This should match N_SIM in main.c

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

    axes.errorbar(N, mean_dev, yerr=CI_95, marker='o', linestyle='-', capsize=5, label=f"r={r}")
    
axes.set_xlabel("Population size N")
axes.set_ylabel(r"$E[\sup_{t \in [0,T(N)]} |X(t)/N - x(t)|]$")
#axes.set_title("Supremum deviation")
axes.legend()
#axes.grid()

plt.show() 