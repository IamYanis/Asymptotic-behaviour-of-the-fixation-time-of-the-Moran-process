/*
  Absorption-time generator for Euler–Maruyama on [0,1]:

    Model A: dY = lambda * sqrt(Y(1-Y)) dt + sqrt(Y(1-Y)) dW
    Model B: dY =                         sqrt(Y(1-Y)) dW

  Absorption rule (discrete):
    - if y_next <= 0 => absorbed at 0
    - if y_next >= 1 => absorbed at 1
    - otherwise continue

  Output (CSV):
    path_id,model,abs_time,abs_state,steps
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

/* ---------------- RNG: xorshift64* + Box-Muller normal ---------------- */

typedef struct {
    uint64_t s;
    int has_spare;
    double spare;
} rng_t;

static uint64_t xorshift64star(uint64_t *x) {
    uint64_t z = *x;
    z ^= z >> 12;
    z ^= z << 25;
    z ^= z >> 27;
    *x = z;
    return z * 2685821657736338717ULL;
}

static void rng_seed(rng_t *r, uint64_t seed) {
    if (seed == 0) seed = 0x9e3779b97f4a7c15ULL;
    r->s = seed;
    r->has_spare = 0;
    r->spare = 0.0;
}

static double rng_u01(rng_t *r) {
    /* uniform in (0,1), avoid exact 0 */
    uint64_t v = xorshift64star(&r->s);
    /* 53-bit mantissa */
    const uint64_t mant = (v >> 11) | 1ULL; /* ensure nonzero */
    return (double)mant * (1.0 / 9007199254740992.0); /* 2^53 */
}

static double rng_n01(rng_t *r) {
    /* Standard normal via Box–Muller (polar form) */
    if (r->has_spare) {
        r->has_spare = 0;
        return r->spare;
    }
    double u, v, s;
    do {
        u = 2.0 * rng_u01(r) - 1.0;
        v = 2.0 * rng_u01(r) - 1.0;
        s = u*u + v*v;
    } while (s >= 1.0 || s == 0.0);
    double mul = sqrt(-2.0 * log(s) / s);
    r->spare = v * mul;
    r->has_spare = 1;
    return u * mul;
}

/* ---------------- SDE simulation ---------------- */

typedef enum {
    MODEL_A = 0, /* with drift lambda */
    MODEL_B = 1  /* no drift */
} model_t;

typedef struct {
    double abs_time;
    int abs_state;   /* 0 or 1 */
    long steps;
    int absorbed;    /* 0/1 */
} absorption_result_t;

static inline double clamp01(double y) {
    if (y < 0.0) return 0.0;
    if (y > 1.0) return 1.0;
    return y;
}

static absorption_result_t simulate_absorption(
    rng_t *rng,
    model_t model,
    double y0,
    double lambda,
    double dt,
    long max_steps
) {
    absorption_result_t res;
    res.abs_time = NAN;
    res.abs_state = -1;
    res.steps = 0;
    res.absorbed = 0;

    double y = y0;

    /* If starting exactly at boundary, absorb immediately */
    if (y <= 0.0) {
        res.abs_time = 0.0;
        res.abs_state = 0;
        res.steps = 0;
        res.absorbed = 1;
        return res;
    }
    if (y >= 1.0) {
        res.abs_time = 0.0;
        res.abs_state = 1;
        res.steps = 0;
        res.absorbed = 1;
        return res;
    }

    const double sqrt_dt = sqrt(dt);

    for (long k = 0; k < max_steps; ++k) {
        double a = y * (1.0 - y);
        if (a < 0.0) a = 0.0; /* numerical guard */
        double sigma = sqrt(a);

        double drift = 0.0;
        if (model == MODEL_A) {
            drift = lambda * sigma;
        }

        double z = rng_n01(rng);
        double y_next = y + drift * dt + sigma * sqrt_dt * z;

        /* Absorption check on crossing */
        if (y_next <= 0.0) {
            res.abs_time = (double)(k + 1) * dt;
            res.abs_state = 0;
            res.steps = k + 1;
            res.absorbed = 1;
            return res;
        }
        if (y_next >= 1.0) {
            res.abs_time = (double)(k + 1) * dt;
            res.abs_state = 1;
            res.steps = k + 1;
            res.absorbed = 1;
            return res;
        }

        /* Keep inside (0,1) if desired; not strictly needed after crossing checks */
        y = clamp01(y_next);
    }

    /* Not absorbed within max_steps */
    res.abs_time = (double)max_steps * dt;
    res.abs_state = -1;
    res.steps = max_steps;
    res.absorbed = 0;
    return res;
}

/* ---------------- CLI ---------------- */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s --npaths N --dt DT --tmax TMAX --y0 Y0 --lambda L --seed SEED\n"
        "\n"
        "Defaults:\n"
        "  N=10000, dt=1e-4, tmax=10, y0=0.5, lambda=0, seed=1\n"
        "\n"
        "Output: CSV to stdout: path_id,model,abs_time,abs_state,steps\n"
        "  model: A (with drift), B (no drift)\n"
        "  abs_state: 0 or 1, or -1 if not absorbed by tmax\n",
        prog
    );
}

static int arg_eq(const char *a, const char *b) { return strcmp(a,b) == 0; }

int main(int argc, char **argv) {
    long npaths = 10000;
    double dt = 1e-4;
    double tmax = 10.0;
    double y0 = 0.5;
    double lambda = 0.0;
    uint64_t seed = 1;

    for (int i = 1; i < argc; ++i) {
        if (arg_eq(argv[i], "--npaths") && i + 1 < argc) npaths = atol(argv[++i]);
        else if (arg_eq(argv[i], "--dt") && i + 1 < argc) dt = atof(argv[++i]);
        else if (arg_eq(argv[i], "--tmax") && i + 1 < argc) tmax = atof(argv[++i]);
        else if (arg_eq(argv[i], "--y0") && i + 1 < argc) y0 = atof(argv[++i]);
        else if (arg_eq(argv[i], "--lambda") && i + 1 < argc) lambda = atof(argv[++i]);
        else if (arg_eq(argv[i], "--seed") && i + 1 < argc) seed = (uint64_t)strtoull(argv[++i], NULL, 10);
        else if (arg_eq(argv[i], "--help") || arg_eq(argv[i], "-h")) { usage(argv[0]); return 0; }
        else {
            fprintf(stderr, "Unknown/invalid arg: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (dt <= 0.0 || tmax <= 0.0 || npaths <= 0) {
        fprintf(stderr, "Invalid parameters.\n");
        return 1;
    }
    if (y0 < 0.0 || y0 > 1.0) {
        fprintf(stderr, "y0 must be in [0,1].\n");
        return 1;
    }

    long max_steps = (long)floor(tmax / dt + 0.5);
    if (max_steps < 1) max_steps = 1;

    rng_t rng;
    rng_seed(&rng, seed);

    printf("path_id,model,abs_time,abs_state,steps\n");

    for (long p = 0; p < npaths; ++p) {
        /* Model A */
        absorption_result_t ra = simulate_absorption(&rng, MODEL_A, y0, lambda, dt, max_steps);
        printf("%ld,A,%.17g,%d,%ld\n", p, ra.abs_time, ra.abs_state, ra.steps);

        /* Model B */
        absorption_result_t rb = simulate_absorption(&rng, MODEL_B, y0, lambda, dt, max_steps);
        printf("%ld,B,%.17g,%d,%ld\n", p, rb.abs_time, rb.abs_state, rb.steps);
    }

    return 0;
}