#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "./main.h"

void print_grid(Grid* g) {
    for (unsigned int i = 0; i < g->width * g->height; i++) {
        printf("%f", g->val[i]);

        if ((i + 1) % g->width != 0) {
            printf(" ");
        } else {
            printf("\n");
        }
    }
}

void generate_grid(Grid* g, int init) {
    unsigned int width = g->width;
    unsigned int height = g->height;
    // Allocate memory for each grid cell
    g->val = (double*)malloc(width * height * sizeof(double));
    if (g->val == NULL) {
        printf("Cannot allocate memory for new grid of size %d, %d", width,
               height);
        exit(-1);
        return;
    }
    unsigned int x, y;
    for (unsigned int i = 0; i < width * height; i++) {
        x = i % width;
        y = i / width;

        // Initialise each memory location with a value
        // If init == True assign left and top edges to 1 and the other
        // edges to 0. The rest are set to random values
        if (init == PATTERN_RANDOM) {
            g->val[i] = (double)(rand() / 10) / (double)(RAND_MAX / 10);
        } else if (init == PATTERN_GRADIENT) {
            if (x == 0 || (y == 0 && x != width - 1)) {
                g->val[i] = 1.0;
            } else if (x == width - 1 || (x != 0 && y == height - 1)) {
                g->val[i] = -1.0;
            } else {
                g->val[i] = 0.0;
            }
        } else if (init == PATTERN_ZERO) {
            g->val[i] = 0;
        }
    }
}

double compute_step(Grid* in, Grid* out) {
    double v0, v1, v2, v3, max_diff = 0.0, diff;
    unsigned int index;
    for (unsigned int y = 1; y < in->height - 1; y++) {
        for (unsigned int x = 1; x < in->width - 1; x++) {
            index = x + y * in->width;
            v0 = in->val[x + 1 + y * in->width];
            v1 = in->val[x - 1 + y * in->width];
            v2 = in->val[x + (y + 1) * in->width];
            v3 = in->val[x + (y - 1) * in->width];
            out->val[index] = (v0 + v1 + v2 + v3) / 4.0;
            diff = fabs(out->val[index] - in->val[index]);
            max_diff = max_diff < diff ? diff : max_diff;
        }
    }
    return max_diff;
}

int main(int argc, char** argv) {
    unsigned int width = 64, height = 64;
    double precision = 0.001;
    if (argc == 4) {
        width = atoi(argv[1]);
        height = atoi(argv[2]);
        precision = atof(argv[3]);
    }

    int is_less_than_prec;
    unsigned int iteration = 0;

    Grid g0, g1, temp;
    g0.width = g1.width = width;
    g0.height = g1.height = height;
    generate_grid(&g0, PATTERN_GRADIENT);
    generate_grid(&g1, PATTERN_GRADIENT);
    while (TRUE) {
        is_less_than_prec = compute_step(&g0, &g1) < precision;
        iteration++;
        if (is_less_than_prec == TRUE) {
            break;
        }
        temp = g1;
        g1 = g0;
        g0 = temp;
    }
    printf("%d,%u,%u,%f,%d,Na,", 1, width, height, precision, iteration);
    print_grid(&g1);
    return 0;
}