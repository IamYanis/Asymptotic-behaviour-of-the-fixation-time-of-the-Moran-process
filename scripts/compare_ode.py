import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ODE function
def odefun(t, y, r):
    dy = ((r - 1) * y * (1 - y)) / (r * y + (1. - y))
    return dy

# Solve ODE using SciPy
def solve_ode_python(T_N, r, x0=0.5):
    t_eval = np.arange(0, T_N, 0.001)  # Time points for evaluation
    odeout = solve_ivp(odefun, [0., T_N], [x0], args=(r,), t_eval=t_eval, method="RK45")
    return odeout.t, odeout.y[0]

# Load data from main.c
file_path = "data/moran_simulation_results.csv"
df = pd.read_csv(file_path)

# Select specific (N, r) for visualization
N_plot = 300
r_plot = 0.8

# Filter ODE solution from main.c
df_ode = df[(df["N"] == N_plot) & (df["r"] == r_plot) & (df["Type"].str.strip() == "ODE")].copy()

if df_ode.empty:
    print(f"Error: No ODE data found for N={N_plot}, r={r_plot}. Check your main.c output.")
    exit()

# Solve ODE using Python's solve_ivp
T_N = np.log(N_plot)
ode_t_python, ode_x_python = solve_ode_python(T_N, r_plot)

# Plot comparison
plt.figure(figsize=(8, 5))
plt.plot(df_ode["Time"], df_ode["Value"], label="ODE (main.c)", linestyle="dashed", color="red")
plt.plot(ode_t_python, ode_x_python, label="ODE (Python, solve_ivp)", linestyle="solid", color="blue")

# Labels and Title
plt.xlabel("Time")
plt.ylabel("Normalized Population (x)")
plt.title(f"Comparison of ODE Solutions (N={N_plot}, r={r_plot})")
plt.legend()
plt.grid()
plt.show()
