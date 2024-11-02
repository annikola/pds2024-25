#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

#define MAX_THREADS 5
#define MAX_CLUSTERS 2

void write_2D_array_to_matfile(const char *filename);
double random_double(double min, double max);
void *distance_calculator(void *ccN);

typedef struct knn_centroid {
    int center_index;
    double **neighbors;
} knn_centroid;

int main(int argc, char *argv[]) {

    int i, j, k, l, t;
    int c_size, q_size, d, closest_centroid;
    double **c, **q;
    double **temp[MAX_CLUSTERS];
    double diff_squares_sum, min_euclid_distance;
    double curr_centroid_center_distances[MAX_CLUSTERS];
    knn_centroid clusters[MAX_CLUSTERS];
    pthread_t thread_ids[MAX_THREADS];

    /* THA FYGEI STO TELOS!!! --> */
    
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

    write_2D_array_to_matfile("test.mat");

    for (i = 0; i < q_size; i++) {
        for (j = 0; j < d; j++) {
            q[i][j] = random_double(1.0, 100.0);
        }
    }

    /* <-- THA FYGEI STO TELOS!!! */

    for (l = 0;  l < MAX_CLUSTERS; l++) {
        clusters[l].neighbors = (double **)malloc(c_size * sizeof(double *));
    }

    for (i = 0; i < MAX_CLUSTERS; i++) {
        clusters[i].center_index = rand() % c_size;
    }

    for (k = 0; k < MAX_CLUSTERS; k++) {
        temp[k] = clusters[k].neighbors;
    }
    for (i = 0; i < c_size; i++) {
        for (k = 0; k < MAX_CLUSTERS; k++) {
            diff_squares_sum = 0.0;
            for (j = 0; j < d; j++) {
                diff_squares_sum += pow(c[i][j] - c[clusters[k].center_index][j], 2);
            }
            curr_centroid_center_distances[k] = sqrt(diff_squares_sum);
        }

        min_euclid_distance = 1e10; // MAX_DOUBLE...change that!
        for (k = 0; k < MAX_CLUSTERS; k++) {
            if (curr_centroid_center_distances[k] < min_euclid_distance) {
                min_euclid_distance = curr_centroid_center_distances[k];
                closest_centroid = k;
            }
        }

        *(temp[closest_centroid]) = c[i];
        temp[closest_centroid]++;
    }

    for (k = 0; k < MAX_CLUSTERS; k++) {
        clusters[k].neighbors = realloc(clusters[k].neighbors, (temp[k] - clusters[k].neighbors) * sizeof(double *));
    }

    // for (k = 0; k < MAX_CLUSTERS; k++) {
    //     for (i = 0; i < temp[k] - clusters[k].neighbors; i++) {
    //         for (j = 0; j < d; j++) {
    //             printf("%lf ", clusters[k].neighbors[i][j]);
    //         }
    //         printf("\b\n");
    //     }
    //     printf("\n");
    // }

    // for (t = 0; t < MAX_THREADS; t++) {
    //     pthread_create(&thread_ids[t], NULL, distance_calculator, &cc2);
    // }

    // for (t = 0; t < MAX_THREADS; t++) {
    //     pthread_join(thread_ids[t], NULL);
    // }

    return 0;
}

void write_2D_array_to_matfile(const char *filename) {

    MATFile *matFile = matOpen(filename, "w");
    if (matFile == NULL) {
        printf("Error opening MAT-file %s\n", filename);
        return;
    }

    double array[3][3] = { {1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0},
                           {7.0, 8.0, 9.0} };

    mxArray *matArray = mxCreateDoubleMatrix(3, 3, mxREAL);
    if (matArray == NULL) {
        printf("Could not create mxArray.\n");
        matClose(matFile);
        return;
    }

    memcpy(mxGetPr(matArray), array, sizeof(array));

    matPutVariable(matFile, "myArray", matArray);

    mxDestroyArray(matArray);
    matClose(matFile);
    printf("2D array written to %s successfully.\n", filename);
}

double random_double(double min, double max) {

    double scale;

    scale = rand() / (double) RAND_MAX;

    return min + scale * (max - min);
}

void *distance_calculator(void *ccN) {
    printf("%d\n", *(int *)ccN);
}
