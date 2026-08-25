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

Numerical simulations for studying the **asymptotic behaviour of the fixation/absorption time of the Moran process** and its deterministic approximation.

This project simulates a continuous-time Moran process for different population sizes and fitness parameters. It compares stochastic trajectories with the corresponding deterministic ODE and investigates quantities such as:

- fixation / absorption times,
- the distribution of absorption times,
- the deviation between the stochastic process and its deterministic approximation,
- the dependence of these quantities on the population size \(N\).

---

## The Moran Process

Consider a population of fixed size \(N\), with \(X_t\) individuals of one type and \(N-X_t\) individuals of the other type.

The process evolves as a continuous-time birth-death process with transition rates

\[q_{X,X+1}=
\frac{(N-X)rX}{rX+(N-X)},
\]

and

\[q_{X,X-1}=\frac{(N-X)X}{rX+(N-X)},\]

where \(r\) is the relative fitness parameter.

The states

\[
X=0
\qquad\text{and}\qquad
X=N
\]

are absorbing.

The normalized population is

\[
x_t = \frac{X_t}{N}.
\]

For large \(N\), the stochastic process can be compared with the deterministic solution of

\[\frac{dx}{dt}=\frac{(r-1)x(1-x)}{rx+(1-x)}.\]

The goal of this repository is to numerically investigate this approximation and the asymptotic behaviour of the process near absorption.

---

## Repository Structure

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
├── Cullen_Frey_plot/
├── plot_Cullen_Frey/
├── plots/
├── include/
├── bin/
│
└── README.md
```

### `src/main.c`

Main simulation program.

It:

- simulates the continuous-time Moran process using exponential waiting times,
- solves the deterministic ODE using the GNU Scientific Library (GSL),
- records stochastic trajectories,
- records absorption times,
- interpolates trajectories on the ODE time grid,
- computes the supremum deviation

\[
\sup_{t\in[0,T_N]}
\left|
\frac{X_t}{N}-x(t)
\right|.
\]

### `scripts/plot_traj.py`

Plots simulated Moran trajectories together with the deterministic ODE solution.

### `scripts/plot_deviation.py`

Studies the mean supremum deviation between the normalized Moran process and its deterministic approximation as a function of \(N\):

\[
\mathbb{E}
\left[
\sup_{t\in[0,T(N)]}
\left|
\frac{X_t}{N}-x(t)
\right|
\right].
\]

The script also provides a logarithmic representation that can be used to investigate the scaling with respect to \(N\).

### `scripts/plot_hist.py`

Plots empirical histograms of the absorption/fixation times for different population sizes.

---

## Requirements

### C

The simulation requires:

- a C compiler supporting C99 or later,
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/).

On Ubuntu/Debian:

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

Install them with:

```bash
pip install numpy pandas matplotlib scipy
```

---

## Compilation

From the root of the repository:

```bash
gcc -O2 -o moran_sim src/main.c -lgsl -lgslcblas -lm
```

On Windows with MSYS2/MinGW, you may need to explicitly specify the GSL include and library directories:

```bash
gcc -O2 -o moran_sim.exe src/main.c \
    -I"C:/msys64/mingw64/include" \
    -L"C:/msys64/mingw64/lib" \
    -lgsl -lgslcblas -lm
```

---

## Running the Simulations

Run:

```bash
./moran_sim
```

or, on Windows:

```bash
./moran_sim.exe
```

Simulation parameters can be configured directly in `src/main.c`.

For example:

```c
#define N_SIM 10000

int N_values[] = {50, 100, 300, 1000};
double r_values[] = {0.8};

double T_N = 500;
int K = 10000;
```

The main parameters are:

| Parameter | Description |
|---|---|
| `N_SIM` | Number of Monte Carlo simulations |
| `N_values` | Population sizes to simulate |
| `r_values` | Relative fitness parameters |
| `T_N` | Simulation time horizon |
| `K` | Number of time subdivisions used for the comparison with the ODE |

---

## Studying the Supremum Deviation

To study

\[
\sup_{t\in[0,T(N)]}
\left|
\frac{X_t}{N}-x(t)
\right|,
\]

choose `K` sufficiently large in `main.c` to obtain a fine discretization of the interval \([0,T(N)]\).

Generate the data:

```bash
./moran_sim
```

Then plot the deviation:

```bash
python scripts/plot_deviation.py
```

Make sure that the value

```python
n_sim = ...
```

in `plot_deviation.py` matches

```c
#define N_SIM ...
```

in `main.c`.

---

## Studying the Fixation-Time Distribution

When only the fixation/absorption-time distribution is needed, the computation of a finely discretized trajectory is unnecessary.

Set:

```c
K = 1;
```

and choose `T_N` large enough to cover the range of fixation times being studied.

For example, for sufficiently large fitness parameters, one may use:

```c
T_N = 500;
```

Generate the simulations:

```bash
./moran_sim
```

and plot the histograms:

```bash
python scripts/plot_hist.py
```

---

## Plotting Moran Trajectories

To compare stochastic Moran trajectories with the deterministic ODE solution:

```bash
python scripts/plot_traj.py
```

This allows the finite-population stochastic dynamics to be visually compared with their deterministic approximation.

---

## Typical Workflow

### 1. Compile the simulation

```bash
gcc -O2 -o moran_sim src/main.c -lgsl -lgslcblas -lm
```

### 2. Run the Monte Carlo simulations

```bash
./moran_sim
```

### 3. Plot stochastic trajectories and the ODE

```bash
python scripts/plot_traj.py
```

### 4. Study the supremum deviation

```bash
python scripts/plot_deviation.py
```

### 5. Study fixation-time distributions

```bash
python scripts/plot_hist.py
```

---

## Numerical Method

The continuous-time Moran process is simulated using exponentially distributed waiting times.

If the total transition rate at the current state is \(\lambda\), the waiting time until the next event is

\[
\Delta t
=
-\frac{\log U}{\lambda},
\qquad
U\sim\operatorname{Uniform}(0,1).
\]

The deterministic approximation is solved numerically using the ODE routines provided by GSL.

---

## Research Objective

The broader goal of this project is to investigate the **large-population asymptotics of the Moran process**, with particular emphasis on its fixation time.

The simulations can be used to study:

1. convergence of the normalized Moran process toward its deterministic approximation,
2. the scaling of the stochastic deviation with \(N\),
3. the distribution of the absorption/fixation time,
4. possible limiting distributions after appropriate centering and scaling.

These numerical experiments provide evidence for understanding the asymptotic behaviour of the model as

\[
N \longrightarrow \infty.
\]

---

## Reproducibility

For reproducible numerical experiments, the pseudo-random number generator can be initialized with a fixed seed.

Using the same seed and simulation parameters allows experiments to be reproduced across runs.

---

## License

No license is currently specified.

If this repository is intended for public reuse or distribution, consider adding a `LICENSE` file.
