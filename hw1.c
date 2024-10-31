#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define MAX_THREADS 5
#define MAX_CLUSTERS 2

double random_double(double min, double max);
void *distance_calculator(void *ccN);

typedef struct knn_centroid {
    int center_index;
    double **neighbors;
} knn_centroid;

int main(int argc, char *argv[]) {

    int i, j, k, t, cc1, cc2;
    int c_size, q_size, d;
    double **c, **q;
    double diff_squares_sum, min_euclid_distance, closest_centroid;
    double curr_centroid_center_distances[MAX_CLUSTERS];
    knn_centroid clusters[MAX_CLUSTERS];
    pthread_t thread_ids[MAX_THREADS];

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

    for (i = 0; i < MAX_CLUSTERS; i++) {
        clusters[i].center_index = rand() % c_size;
    }

    cc1 = cc2 = 0;
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
        // clusters->neighbors Ylopoiise to...
        if (closest_centroid == 0){
            cc1++;
        }
        else {
            cc2++;
        }
    }
    printf("cc1: %d\ncc2: %d\n", cc1, cc2);

    for (t = 0; t < MAX_THREADS; t++) {
        pthread_create(&thread_ids[t], NULL, distance_calculator, &cc2);
    }

    for (t = 0; t < MAX_THREADS; t++) {
        pthread_join(thread_ids[t], NULL);
    }

    // for (i = 0; i < c_size; i++) {
    //     for (j = 0; j < d; j++) {
    //         printf("%lf ", c[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n");

    // for (i = 0; i < q_size; i++) {
    //     for (j = 0; j < d; j++) {
    //         printf("%lf ", q[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n");

    return 0;
}

double random_double(double min, double max) {

    double scale;

    scale = rand() / (double) RAND_MAX;

    return min + scale * (max - min);
}

void *distance_calculator(void *ccN) {
    printf("%d\n", *(int *)ccN);
}
