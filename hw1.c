#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double random_double(double min, double max);

int main(int argc, char *argv[]) {

    int i, j;
    int c_size, q_size, d;
    double **c, **q;

    if (argc < 4) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    srand(time(0));

    c_size = atoi(argv[1]);
    q_size = atoi(argv[2]);
    d = atoi(argv[3]);

    c = (double **)malloc(c_size * sizeof(double *));
    for (i = 0; i < c_size; i++) {
        c[i] = (double *)malloc(d * sizeof(double));
    }

    q = (double **)malloc(q_size * sizeof(double *));
    for (j = 0; j < q_size; j++) {
        q[j] = (double *)malloc(d * sizeof(double));
    }

    for (i = 0; i < c_size; i++) {
        for (j = 0; j < d; j++) {
            c[i][j] = random_double(1.0, 100.0);
        }
    }

    for (i = 0; i < q_size; i++) {
        for (j = 0; j < d; j++) {
            q[i][j] = random_double(1.0, 100.0);
        }
    }

    for (i = 0; i < c_size; i++) {
        for (j = 0; j < d; j++) {
            printf("%lf ", c[i][j]);
        }
        printf("\b\n");
    }
    printf("\n");

    for (i = 0; i < q_size; i++) {
        for (j = 0; j < d; j++) {
            printf("%lf ", q[i][j]);
        }
        printf("\b\n");
    }
    printf("\n");

    return 0;
}

double random_double(double min, double max) {

    double scale;

    scale = rand() / (double) RAND_MAX;

    return min + scale * (max - min);
}
