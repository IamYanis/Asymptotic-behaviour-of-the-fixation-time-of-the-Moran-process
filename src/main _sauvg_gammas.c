#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PROGRESS_FILE "data/progress_state.txt"

#ifdef _WIN32
    #include <direct.h>  // For _mkdir() on Windows
#else
    #include <sys/stat.h>  // For mkdir() on Linux/Mac
#endif

//#define N_SIM 10000   // Number of MC simulations
#define N_SIM 10000
//#define DT 0.001       // Time step for ODE integration
//int K = (int)(10000);  //scalling precision in the set [0, T(N)] (FOR COMPUTING THE SUPREMUM DEVIATION)
//int K = (int)(1); //(USELESS IF LOOKING AT THE TIME OF ABSORPTION)

// Function to compute birth rate
double rplus(int X, int N, double r) {
    return ((N - X) * r * X) / (r * X + (N - X));
}

// Function to compute death rate
double rminus(int X, int N, double r) {
    return ((N - X) * X) / (r * X + (N - X));
}

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

// Define the logistic ODE function
int logistic_ode(double t, const double y[], double dydt[], void *params) {
    double r = *(double *)params;
    dydt[0] = ((r - 1) * y[0] * (1 - y[0])) / (r * y[0] + (1 - y[0]));
    return GSL_SUCCESS;
}

// Solve logistic ODE using GSL
void solve_ode_gsl(double *ode_t, double *ode_x, int steps, double r, double x0, double DT) {
    gsl_odeiv2_system sys = {logistic_ode, NULL, 1, &r};

    // Adaptive Runge-Kutta-Fehlberg (RKF45) solver
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, DT, 1e-6, 1e-6);

    double t = 0, y[1] = {x0};
    
    for (int i = 0; i < steps; i++) {
        double ti = i * DT;
        gsl_odeiv2_driver_apply(driver, &t, ti, y);
        ode_t[i] = t;
        ode_x[i] = y[0];
    }

    gsl_odeiv2_driver_free(driver);
}


// Interpolate Moran process at a given ODE time
double interpolate_moran(double *moran_t, double *moran_x, int count, double time) {
    if (time <= moran_t[0]) return moran_x[0];
    if (time >= moran_t[count - 1]) return moran_x[count - 1];

    for (int i = 0; i < count - 1; i++) {
        if (moran_t[i] <= time && time <= moran_t[i + 1]) {
            double denominator = moran_t[i + 1] - moran_t[i];
            if (denominator == 0) return moran_x[i];  // Handle division by zero

            double weight = (time - moran_t[i]) / denominator;

            // Debugging print statement
            if (weight == 0.0) {
                //printf("Weight is 0 at index %d, time: %f, moran_t[%d]: %f\n", i, time, i, moran_t[i]);
            } else if (weight == 1.0) {
                //printf("Weight is 1 at index %d, time: %f, moran_t[%d+1]: %f\n", i, time, i, moran_t[i + 1]);
            }

            return (1 - weight) * moran_x[i] + weight * moran_x[i + 1];
        }
    }

    return moran_x[count - 1]; // Return last value if out of bounds
}


// Function to simulate Moran process
void simulate_moran_process() {
    // Ensure directory exists
#ifdef _WIN32
    _mkdir("data");  // Windows version
#else
    mkdir("data", 0700);  // Linux/Mac version
#endif


    int next_i = 0, next_j = 0, next_sim = 0;
    int resume = 0;

    // Try to read progress file
    FILE *state_in = fopen(PROGRESS_FILE, "r");
    if (state_in) {
        if (fscanf(state_in, "%d %d %d", &next_i, &next_j, &next_sim) == 3) {
            resume = 1;
            printf("Resuming from i=%d, j=%d, sim=%d\n", next_i, next_j, next_sim);
        }
        fclose(state_in);
    } else {
        printf("No progress file found, starting fresh.\n");
    }


    // Output file
    FILE *abs_out;
    if (resume) {
        // Append if resuming
        abs_out = fopen("data/absorption_times.csv", "a");
    } else {
        // Fresh start: overwrite + header
        abs_out = fopen("data/absorption_times.csv", "w");
        if (abs_out) {
            fprintf(abs_out, "N,gamma,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed\n");
        }
    }

    if (!abs_out) {
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }

    // Write header
    //fprintf(file_out, "N,gamma,r,Type,Sim_ID,Time,Value,std_sup_dev\n");
    //fprintf(abs_out, "N,gamma,r,Sim_ID,AbsorptionTime,AbsorbingState,Absorbed\n");

    int N_values[] = {1000};
    //int N_values[] = {10, 20, 30, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300};
    //int N_values[] = {500};
    //int N_values[] = {10, 20, 50, 70, 100, 150, 200, 250, 300};//, 350, 400, 450, 500};//, 550, 600, 650, 700, 750, 800, 850, 900, 950,1000}; //, 1250, 1500};
    //int N_values[] = {10, 20, 50, 70, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900,1000};
    //int N_values[] = {1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    //int N_values[] = {10, 20, 50, 70, 100, 150, 200, 250, 300, 350, 400, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500};
    //int N_values[] = {500, 700, 1000, 1500, 2000, 3000};
    //double r_values[] = {1.05, 1.1, 1.2, 1.5};


//double gamma_values[] = {.4, .5, .7, .8, .9, 1.0, 1.1, 1.2};
double gamma_values[] = {.5, .6, .7, .8, .9, 1.0, 1.1, 1.2};
//double gamma_values[] = {.5,.525, .55,.575, .6,.625, .65,.675, .7, .725, .75,.775, .8,.825, .85,.875,  .9, .925, .95,.975, 1.0, 1.025, 1.05,1.075, 1.1,1.125, 1.15, 1.175, 1.2,1.225, 1.25,1.275, 1.3,1.325, 1.35,1.375, 1.4,1.425, 1.45,1.475, 1.5};

    int N_size = sizeof(N_values)/sizeof(N_values[0]);
    int r_size = sizeof(gamma_values)/sizeof(gamma_values[0]); //number of gammas

    //double r_values[] = {0.8};

    //double r_values[] = {0.9};
    //double r_values[] = {1.3, 1.4};
    //double r_values[] = {1.01, 1.03, 1.05, 1.07, 1.1, 1.2, 1.3, 1.5};
    //double r_values[] = {1.3, 1.2, 1.1, 1.03};
    //int N_size = sizeof(N_values) / sizeof(N_values[0]);
    //int r_size = sizeof(r_values) / sizeof(r_values[0]);


    srand(645321);
    //srand(time(NULL));
    
    //double hh = 0.3;

    for (int i = next_i; i < N_size; i++) {
        int j_start = (i == next_i ? next_j : 0);

        for (int j = j_start; j < r_size; j++) {
            int sim_start = (i == next_i && j == next_j ? next_sim : 0);

            int N = N_values[i];
            double gamma = gamma_values[j];
            //double r = r_values[j];
            double r = 1.0 - pow(N, -gamma);
            double T_N = 10000;
            //double T_N = 100;
            //double T_N = log(N);
            double res_fitness = fabs(r-1);
            //double T_N = 1/(res_fitness)*floor(exp( (log(2)*log(2)*pow(N, hh))/(2*res_fitness) ));
            printf("|1-r|=%.6f\n", res_fitness);
            printf("T_N=%.6f\n", T_N);
            printf("gamma=%.6f\n",gamma);

            int K = (int)(10000);  //scalling precision in the set [0, T(N)] (FOR COMPUTING THE SUPREMUM DEVIATION)
            //int K = (int)(1); //(USELESS IF LOOKING AT THE TIME OF ABSORPTION)
            double DT_adaptive = T_N / K;
            int steps = (int)(T_N / DT_adaptive) + 1;

            //printf("Processing N=%d, gamma=%.3f, r=%.2f...\n", N, gamma, r);
            printf("Processing N=%d, r=%.2f...\n", N, r);
            printf("|r-1|=%.6f\n", fabs(r-1.0));
            printf("T_N=%.6f\n", T_N);

            // Allocate memory for ODE solution
            double *ode_t = (double *)malloc(steps * sizeof(double));
            double *ode_x = (double *)malloc(steps * sizeof(double));
            if (!ode_t || !ode_x) {
                perror("Memory allocation failed");
                exit(EXIT_FAILURE);
            }

            // Compute and store ODE solution
            //printf("  Computing ODE solution...\n");
            // Compute and store ODE solution using GSL
            solve_ode_gsl(ode_t, ode_x, steps, r, 0.5, DT_adaptive);

            // Don't need to store 
            //for (int k = 0; k < steps; k++) {
            //    fprintf(file_out, "%d,%.6f,%.6f,ODE,-1,%.10f,%.10f,\n", N, gamma, r, ode_t[k], ode_x[k]);
            //}

            //printf("  ODE solution stored.\n");

            // Simulate Moran process
            //printf("  Running Moran simulations (%d simulations)...\n", N_SIM);
            double sup_deviations[N_SIM];
            double total_sup_deviation = 0.0;
            double total_sup_deviation_sq = 0.0;

            for (int sim = sim_start; sim < N_SIM; sim++) {
                //int X = N / 2;
                int X = N/2;
                //int X = N-1; 
                //int X = N/2;
                //int X = .7 * N; 
                //int X = (int)(N - sqrt(N));
                //int X = (int)(N - pow(N, 0.25));
                double t = 0.0;
                double *moran_t = (double *)malloc(steps * sizeof(double));
                double *moran_x = (double *)malloc(steps * sizeof(double));
                int moran_count = 0;
                double max_dev = 0.0;

                double last_valid_t = t;  // Keep track of last valid time

                while (t < T_N && X > 0 && X < N) {
                    double birth_rate = rplus(X, N, r);
                    double death_rate = rminus(X, N, r);
                    double total_rate = birth_rate + death_rate;

                    if (total_rate < 1e-6) {  // Threshold to detect small rates
                        printf("Warning: Small total_rate detected! X = %d, t = %.6f, total_rate = %.12f\n", X, t, total_rate);
                    }

                    if (total_rate <= 0) break; // Prevent infinite loop


                    // Initial computation of rand_val but bugged because often be zero:
                    //double rand_val = (double)rand() / RAND_MAX;

                    // Alternative way that prevent from being zero but can introduce a bias in the law
                    //double rand_val;
                    //do {
                    //    rand_val = (double)rand() / RAND_MAX;  // Generate number in (0,1]
                    //} while (rand_val == 0);  // Reject 0 and retry

                    // Alternative way that impose to rand_val a minimum value: 
                    double rand_val = (double)rand() / RAND_MAX;

                    // Impose a minimal value if rand_val is too small
                    const double MIN_RAND = 1e-12;  // Smallest allowed value
                    if (rand_val < MIN_RAND) {
                        rand_val = MIN_RAND;
}

                    double log_val = -log(rand_val) / total_rate;

                    if (rand_val == 0) {
                        printf("Warning: rand() produced 0! This will cause Inf in log().\n");
                    }
                
                    // Check if log calculation goes bad
                    if (!isfinite(log_val)) {
                        printf("Warning: log() calculation resulted in Inf or NaN! rand_val = %.12f, log_val = %.12f\n",
                               rand_val, log_val);
                    }
                
                    t += log_val;  // Update time

                    //t += -log((double)rand() / RAND_MAX) / total_rate;

                    if (!isfinite(t)) {  // Check for Inf or NaN
                        //printf("Warning: Inf absorption time, using last valid time t=%.6f\n", last_valid_t);
                        t = last_valid_t;  // Use the last valid time instead of Inf
                        break;
                    }
                
                    last_valid_t = t;  // Update last valid time

                    if ((double)rand() / RAND_MAX < (birth_rate / total_rate)) {
                        X += 1;
                    } else {
                        X -= 1;
                    }

                    // Store Moran process trajectory
                    if (moran_count < steps) {
                        moran_t[moran_count] = t;
                        moran_x[moran_count] = X / (double)N;
                        moran_count++;
                    }
                    
                    // Don't need to be stored
                    //fprintf(file_out, "%d,%.6f,%.6f,Moran,%d,%.10f,%.10f,\n", N, gamma, r, sim, t, X / (double)N);

                }

                // progress update
                printf("\rSimulation %d / %d      ", sim + 1, N_SIM);  // a few spaces wipe leftovers
                fflush(stdout);


                // Determine absorbing state
                int absorbed = (X == 0 || X == N) ? 1 : 0;
                int absorbing_state = absorbed ? ((X == 0) ? 0 : N) : -1; //-1 for not absorbed

                // Save absorption time
                fprintf(abs_out, "%d,%.6f,%.6f,%d,%.10f,%d,%d\n", N, gamma, r, sim, t, absorbing_state, absorbed);
                //fprintf(abs_out, "%d,%.6f,%d,%.10f,%d,%d\n", N, r, sim, t, absorbing_state, absorbed);
                if ((sim % 100) == 0) fflush(abs_out);

  


                
                // Compute supremum deviation over **all** ODE times
                for (int k = 0; k < steps; k++) {
                    double moran_interp = interpolate_moran(moran_t, moran_x, moran_count, ode_t[k]);
                    double deviation = fabs(moran_interp - ode_x[k]);
                    if (deviation > max_dev) {
                        max_dev = deviation;
                    }
                }

            
                // Store supremum deviation for this simulation
                sup_deviations[sim] = max_dev;
                total_sup_deviation += max_dev;
                total_sup_deviation_sq += max_dev * max_dev;

                free(moran_t);
                free(moran_x);

                // === Update progress file to point to the NEXT task ===
                int ni = i, nj = j, ns = sim + 1;

                // If we've finished all sims for this (i, j), move to the next (i, j)
                if (ns >= N_SIM) {
                    ns = 0;
                    nj = j + 1;
                    if (nj >= r_size) {
                        nj = 0;
                        ni = i + 1;
                    }
                }

                // Write next_i, next_j, next_sim to file
                FILE *state_out = fopen(PROGRESS_FILE, "w");
                if (state_out) {
                    fprintf(state_out, "%d %d %d\n", ni, nj, ns);
                    fclose(state_out);
                }
            }

            printf("\n"); // To delete the progression line every time

            // Store mean supremum deviation
            double mean_sup_dev = total_sup_deviation / N_SIM;
            double variance = (total_sup_deviation_sq / N_SIM) - (mean_sup_dev * mean_sup_dev);
            double std_sup_dev = (variance > 0) ? sqrt(variance) : 0.0;

            // Don't need to be stored
            //fprintf(file_out, "%d,%.6f,%.6f,Deviation,-1,%.10f,%.10f,%.10f\n", N, gamma, r, T_N, mean_sup_dev, std_sup_dev);
            //printf("  Mean supremum deviation for N=%d, gamma=%.3f, r=%.2f: %f\n", N, gamma, r, mean_sup_dev);

            // Free memory (AFTER ALL SIMULATIONS)
            free(ode_t);
            free(ode_x);
        }
    }

    //fclose(file_out);
    fclose(abs_out);

    // All done successfully → delete progress file so a new run starts fresh
    remove(PROGRESS_FILE);

    printf("Moran model simulations completed and stored in 'data/absorption_times.csv'.\n");
    }

int main() {
    simulate_moran_process();
    return 0;
}
