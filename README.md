# Moran-project

1. (i) To study the supremum deviation: 
  -> Choose the size of the subdivision K of the set [0,T(N)] large enough

2. (i) To study histograms of the law of the fixation time:
  -> In main.c, choose K=1 to avoid useless computation and T_N big enough to cover the value studied of the law of the fixation time (eg. take T_N = 500 for r >= 1.1)

To create the data: 
  gcc -o moran_sim src/main.c -I"C:/msys64/mingw64/include" -L"C:/msys64/mingw64/lib" -lgsl -lgslcblas -lm

  ./moran_sim

1. (ii) To plot the deviation:
   -> in plot_devation.py, adapt n_sim to match N_SIM in main.c

   python scripts/plot_deviation.py

2. (ii) To plot the histograms of the fixation time:
   python scripts/plot_hist.py

3. To plot the trajectories and ODE:
  python scripts/plot_traj.py


----------------------------------------------------------------------------------

# Asymptotic Behaviour of the Fixation Time of the Moran Process

Numerical simulations for studying the **asymptotic behaviour of the fixation/absorption time of the Moran process**, together with its deterministic approximation.

The project simulates a continuous-time Moran process for different population sizes and fitness parameters, compares stochastic trajectories with the corresponding logistic ODE, and investigates quantities such as:

* fixation / absorption times,
* the distribution of absorption times,
* the deviation between the stochastic process and its deterministic approximation,
* the dependence of these quantities on the population size (N).

## Moran process

Consider a population of fixed size (N), with (X_t) individuals of one type and (N-X_t) individuals of the other type.

The process evolves as a continuous-time birth-death process with transition rates

[
q_{X,X+1}
=========

\frac{(N-X)rX}{rX+(N-X)},
]

and

[
q_{X,X-1}
=========

\frac{(N-X)X}{rX+(N-X)},
]

where (r) is the relative fitness parameter.

The states

[
X=0
\qquad\text{and}\qquad
X=N
]

are absorbing.

The corresponding normalized population is

[
x_t=\frac{X_t}{N}.
]

For large (N), the stochastic process can be compared with the deterministic solution of

[
\frac{dx}{dt}
=============

\frac{(r-1)x(1-x)}
{rx+(1-x)}.
]

This repository numerically explores this approximation and the behaviour of the process near absorption.

## Repository structure

```text
.
├── src/
│   └── main.c
│
├── scripts/
│   ├── plot_deviation.py
│   ├── plot_hist.py
│   └── plot_traj.py
│
├── data/                   # generated simulation data
├── plots/                  # generated figures
├── Cullen_Frey_plot/       # Cullen–Frey analysis
├── plot_Cullen_Frey/
├── include/
├── bin/
│
└── README.md
```

### `src/main.c`

Main simulation program.

It:

* simulates the continuous-time Moran process using exponential waiting times,
* solves the deterministic ODE using the GNU Scientific Library (GSL),
* records stochastic trajectories,
* records absorption times,
* interpolates trajectories on the ODE time grid,
* computes the supremum deviation

[
\sup_{t\in[0,T_N]}
\left|
\frac{X_t}{N}-x(t)
\right|.
]

The generated files are

```text
data/moran_simulation_results.csv
data/absorption_times.csv
```

### `scripts/plot_traj.py`

Plots simulated Moran trajectories together with the deterministic ODE solution.

The population size, fitness parameter, and displayed time interval can be configured directly in the script:

```python
N_plot = 300
r_plot = 0.8
time_range = (0, 30)
```

### `scripts/plot_deviation.py`

Studies the mean supremum deviation between the normalized Moran process and its deterministic approximation as a function of (N).

It produces both the original and logarithmic representations of

[
\mathbb{E}
\left[
\sup_{t\in[0,T(N)]}
\left|
\frac{X_t}{N}-x(t)
\right|
\right].
]

The script also estimates a scaling exponent using linear regression.

### `scripts/plot_hist.py`

Plots histograms of the simulated absorption times for several population sizes.

This can be used to investigate how the distribution of the fixation/absorption time changes as (N) increases.

## Requirements

### C

The simulation requires:

* a C compiler supporting C99 or later,
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/).

On Linux, for example, GSL can typically be installed with

```bash
sudo apt install libgsl-dev
```

On macOS with Homebrew:

```bash
brew install gsl
```

### Python

The plotting scripts require Python 3 and the following packages:

```text
numpy
pandas
matplotlib
scipy
```

Install them with

```bash
pip install numpy pandas matplotlib scipy
```

## Compilation

From the root of the repository, compile the simulation with

```bash
gcc -O2 -o moran_sim src/main.c -lgsl -lgslcblas -lm
```

On Windows with MSYS2/MinGW, you may need to explicitly provide the GSL include and library directories:

```bash
gcc -O2 -o moran_sim.exe src/main.c \
    -I"C:/msys64/mingw64/include" \
    -L"C:/msys64/mingw64/lib" \
    -lgsl -lgslcblas -lm
```

## Running the simulations

Run

```bash
./moran_sim
```

or, on Windows,

```bash
./moran_sim.exe
```

The program automatically creates the `data/` directory if necessary.

Simulation parameters are currently configured directly in `src/main.c`.

For example:

```c
#define N_SIM 10000

int N_values[] = {50, 100, 300, 1000};
double r_values[] = {0.8};

double T_N = 500;
int K = 10000;
```

Here:

* `N_SIM` is the number of Monte Carlo simulations,
* `N_values` contains the population sizes,
* `r_values` contains the fitness parameters,
* `T_N` determines the simulation time horizon,
* `K` determines the discretization used when comparing the stochastic process with the ODE.

## Plotting the results

After running the C simulation, the Python scripts can be executed from the repository root.

### Moran trajectories vs. deterministic limit

```bash
python scripts/plot_traj.py
```

This displays several stochastic trajectories together with the logistic ODE solution.

### Supremum deviation

```bash
python scripts/plot_deviation.py
```

Make sure that

```python
n_sim = ...
```

in `plot_deviation.py` agrees with `N_SIM` in `src/main.c` before computing confidence intervals.

### Absorption-time distributions

```bash
python scripts/plot_hist.py
```

The script reads

```text
data/absorption_times.csv
```

and plots empirical absorption-time histograms for selected values of (N).

## Output data

### `moran_simulation_results.csv`

The main output file has columns of the form

```text
N,r,Type,Sim_ID,Time,Value,std_sup_dev
```

The `Type` column distinguishes between:

* `Moran` — stochastic Moran trajectories,
* `ODE` — deterministic ODE solution,
* `Deviation` — summary statistics for the supremum deviation.

### `absorption_times.csv`

Absorption times are stored as

```text
N,r,Sim_ID,AbsorptionTime
```

and are used for the empirical study of the fixation-time distribution.

## Typical workflow

```bash
# 1. Compile
gcc -O2 -o moran_sim src/main.c -lgsl -lgslcblas -lm

# 2. Run Monte Carlo simulations
./moran_sim

# 3. Compare stochastic trajectories with the ODE
python scripts/plot_traj.py

# 4. Study the supremum deviation
python scripts/plot_deviation.py

# 5. Study absorption-time distributions
python scripts/plot_hist.py
```

## Reproducibility

The C implementation currently initializes the pseudo-random number generator with a fixed seed:

```c
srand(123456);
```

Consequently, simulations are reproducible when using the same implementation and environment.

To generate a different sequence at each execution, this can be replaced by

```c
srand(time(NULL));
```

## Numerical considerations

The simulation uses exponentially distributed waiting times. If the current total transition rate is (\lambda), the next waiting time is generated as

[
\Delta t
========

-\frac{\log U}{\lambda},
\qquad
U\sim\mathrm{Uniform}(0,1).
]

The implementation prevents (U) from becoming numerically too close to zero before evaluating the logarithm.

The deterministic equation is solved using the **Runge–Kutta–Fehlberg RKF45** solver provided by GSL.

## Research objective

The broader goal of the project is to numerically investigate the large-population behaviour of the Moran process and, in particular, the asymptotic behaviour of its fixation time.

The simulations provide a way to compare finite-(N) stochastic dynamics with deterministic approximations and to investigate the scaling and limiting distribution of quantities associated with absorption.

## License

No license is currently specified for this repository. If the project is intended to be publicly reused or distributed, consider adding a `LICENSE` file.
