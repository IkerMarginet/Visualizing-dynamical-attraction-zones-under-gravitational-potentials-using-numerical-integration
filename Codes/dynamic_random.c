#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define MAX_ATTRACTORS 10

typedef struct {
    double k;
    double x, y;
    double color[3];
} Attractor;

typedef struct {
    double x, y;
} Vec2;

// -------------------- Vector Tools --------------------
Vec2 vec_sub(Vec2 a, Vec2 b) { return (Vec2){a.x - b.x, a.y - b.y}; }
Vec2 vec_add(Vec2 a, Vec2 b) { return (Vec2){a.x + b.x, a.y + b.y}; }
Vec2 vec_scale(Vec2 v, double s) { return (Vec2){v.x * s, v.y * s}; }
double vec_norm(Vec2 v) { return sqrt(v.x * v.x + v.y * v.y); }

// -------------------- Gravitational Force --------------------
Vec2 force_on_particle(Vec2 pos, Attractor* attractors, int n) {
    Vec2 f = {0, 0};
    for (int i = 0; i < n; ++i) {
        Vec2 r_vec = vec_sub(pos, (Vec2){attractors[i].x, attractors[i].y});
        double r = vec_norm(r_vec);
        if (r < 1e-9) continue;
        double scale = -attractors[i].k / (r * r * r);
        f = vec_add(f, vec_scale(r_vec, scale));
    }
    return f;
}

// -------------------- RK4 Integrator --------------------
int integrate_rk4(Vec2 pos, Vec2 vel, Attractor* attractors, int n, double dt, int n_steps, double r_stop) {
    for (int step = 0; step < n_steps; ++step) {
        Vec2 a1 = force_on_particle(pos, attractors, n);
        Vec2 k1v = vec_scale(a1, dt);
        Vec2 k1p = vec_scale(vel, dt);

        Vec2 a2 = force_on_particle(vec_add(pos, vec_scale(k1p, 0.5)), attractors, n);
        Vec2 k2v = vec_scale(a2, dt);
        Vec2 k2p = vec_scale(vec_add(vel, vec_scale(k1v, 0.5)), dt);

        Vec2 a3 = force_on_particle(vec_add(pos, vec_scale(k2p, 0.5)), attractors, n);
        Vec2 k3v = vec_scale(a3, dt);
        Vec2 k3p = vec_scale(vec_add(vel, vec_scale(k2v, 0.5)), dt);

        Vec2 a4 = force_on_particle(vec_add(pos, k3p), attractors, n);
        Vec2 k4v = vec_scale(a4, dt);
        Vec2 k4p = vec_scale(vec_add(vel, k3v), dt);

        vel = vec_add(vel, vec_scale(vec_add(vec_add(k1v, vec_scale(k2v, 2)), vec_add(vec_scale(k3v, 2), k4v)), 1.0 / 6.0));
        pos = vec_add(pos, vec_scale(vec_add(vec_add(k1p, vec_scale(k2p, 2)), vec_add(vec_scale(k3p, 2), k4p)), 1.0 / 6.0));

        for (int i = 0; i < n; ++i)
            if (vec_norm(vec_sub(pos, (Vec2){attractors[i].x, attractors[i].y})) < r_stop) return i;
        if (vec_norm(pos) > 2.0) return -1;
    }
    return -1;
}

// -------------------- Symplectic Integrator --------------------
int integrate_symplectic(Vec2 pos, Vec2 vel, Attractor* attractors, int n, double dt, int n_steps, double r_stop) {
    for (int step = 0; step < n_steps; ++step) {
        vel = vec_add(vel, vec_scale(force_on_particle(pos, attractors, n), 0.5 * dt));
        pos = vec_add(pos, vec_scale(vel, dt));
        vel = vec_add(vel, vec_scale(force_on_particle(pos, attractors, n), 0.5 * dt));

        for (int i = 0; i < n; ++i)
            if (vec_norm(vec_sub(pos, (Vec2){attractors[i].x, attractors[i].y})) < r_stop) return i;
        if (vec_norm(pos) > 2.0) return -1;
    }
    return -1;
}

// -------------------- Map Generation --------------------
void generate_map(Attractor* attractors, int n, const char* integrator_name, int grid_size, double dt, int n_steps, double r_stop, const char* filename) {
    unsigned char* image = malloc(3 * grid_size * grid_size);
    Vec2 vel0 = {0, 0};

    int (*integrator)(Vec2, Vec2, Attractor*, int, double, int, double) = NULL;

    if (strcmp(integrator_name, "rk4") == 0)
        integrator = integrate_rk4;
    else if (strcmp(integrator_name, "symplectic") == 0)
        integrator = integrate_symplectic;
    else {
        fprintf(stderr, "Unknown integrator.\n");
        exit(1);
    }

    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            double x = -1.0 + 2.0 * j / (grid_size - 1);
            double y = -1.0 + 2.0 * i / (grid_size - 1);
            Vec2 pos = {x, y};

            int idx = integrator(pos, vel0, attractors, n, dt, n_steps, r_stop);
            unsigned char* pixel = &image[3 * (i * grid_size + j)];

            if (idx >= 0) {
                pixel[0] = (unsigned char)(255 * attractors[idx].color[0]);
                pixel[1] = (unsigned char)(255 * attractors[idx].color[1]);
                pixel[2] = (unsigned char)(255 * attractors[idx].color[2]);
            } else {
                pixel[0] = pixel[1] = pixel[2] = 255;
            }
        }
    }

    FILE* f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", grid_size, grid_size);
    fwrite(image, 1, 3 * grid_size * grid_size, f);
    fclose(f);
    free(image);
}

// -------------------- Random Double --------------------
double rand_double(double min, double max) {
    return min + (max - min) * rand() / (double)RAND_MAX;
}

// -------------------- Main --------------------
int main() {
    srand(time(NULL));

    int n = 2 + rand() % (MAX_ATTRACTORS - 1); // between 2 and MAX_ATTRACTORS
    Attractor attractors[MAX_ATTRACTORS];

    for (int i = 0; i < n; ++i) {
        attractors[i].k = rand_double(0.5, 2.0);
        attractors[i].x = rand_double(-1.0, 1.0);
        attractors[i].y = rand_double(-1.0, 1.0);
        attractors[i].color[0] = rand_double(0.2, 1.0);
        attractors[i].color[1] = rand_double(0.2, 1.0);
        attractors[i].color[2] = rand_double(0.2, 1.0);
    }

    int grid_size = 500;
    double dt = 0.004;
    int n_steps = 5000;
    double r_stop = 0.03;

    generate_map(attractors, n, "rk4", grid_size, dt, n_steps, r_stop, "random_rk4.ppm");
    generate_map(attractors, n, "symplectic", grid_size, dt, n_steps, r_stop, "random_symplectic.ppm");

    printf("Generated images: random_rk4.ppm and random_symplectic.ppm with %d attractors.\n", n);
    return 0;
}
