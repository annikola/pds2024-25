#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <cblas.h>

#define MAX_SET_SPLIT 2
#define Q_SPLIT 4
#define MAX_GB 1

typedef struct {
    double distance;
    int index;
} Neighbor;

typedef struct {
    double *coordinates;
    Neighbor *neighbors;
    int index;
} Point;

typedef struct {
    double *corpus;
    double *query_part;
    Point **query_part_points;
    int corpus_size;
    int query_part_size;
    int dimensions;
    int knns;
} calculate_distances_args;

int qsort_compare(const void* a, const void* b);
double *calculate_norms(double *matrix2D, int m_size, int d);
void *calculate_distances(void *args);
void knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **set_points);

int qsort_compare(const void *a, const void *b) {
    const Neighbor *neighbor1 = (const Neighbor *)a;
    const Neighbor *neighbor2 = (const Neighbor *)b;

    if (neighbor1->distance < neighbor2->distance) {
        return -1; // neighbor1 is less than neighbor2
    } else if (neighbor1->distance > neighbor2->distance) {
        return 1;  // neighbor1 is greater than neighbor2
    } else {
        return 0;  // they are equal
    }
}

double *calculate_norms(double *matrix2D, int m_size, int d) {

    int i, j;
    double *M_norms;

    M_norms = (double *)malloc(m_size * sizeof(double));
    for (i = 0; i < m_size; i++) {
        M_norms[i] = 0;
        for (j = 0; j < d; j++) {
            M_norms[i] += matrix2D[i * d + j] * matrix2D[i * d + j];
        }
    }

    return M_norms;
}

void *calculate_distances(void *args) {

    int i, j, k;
    double *distances_part_matrix, *c_norms, *q_norms, *matrixmul;
    calculate_distances_args *cd_args;

    cd_args = (calculate_distances_args *)args;

    matrixmul = (double *)malloc(cd_args->corpus_size * cd_args->query_part_size * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, cd_args->query_part_size, cd_args->corpus_size, cd_args->dimensions, -2.0, cd_args->query_part, cd_args->dimensions, cd_args->corpus, cd_args->dimensions, 0.0, matrixmul, cd_args->corpus_size);

    c_norms = calculate_norms(cd_args->corpus, cd_args->corpus_size, cd_args->dimensions);
    q_norms = calculate_norms(cd_args->query_part, cd_args->query_part_size, cd_args->dimensions);

    distances_part_matrix = (double *)malloc(cd_args->query_part_size * cd_args->corpus_size * sizeof(double));
    for (i = 0; i < cd_args->query_part_size; i++) {
        for (j = 0; j < cd_args->corpus_size; j++) {
            distances_part_matrix[i * cd_args->corpus_size + j] = q_norms[i] + c_norms[j] + matrixmul[i * cd_args->corpus_size + j];
            
            // Check if the value to be squared is slightly below zero...(UNDERFLOW!!!)
            if (distances_part_matrix[i * cd_args->corpus_size + j] < 0) {
                distances_part_matrix[i * cd_args->corpus_size + j] = 0.0;
            } else {
                distances_part_matrix[i * cd_args->corpus_size + j] = sqrt(distances_part_matrix[i * cd_args->corpus_size + j]);
            }
        }
    }

    for (i = 0; i < cd_args->query_part_size; i++) {
        // if STITCHING REALLOC!!!
        cd_args->query_part_points[i]->neighbors = (Neighbor *)malloc(cd_args->corpus_size * sizeof(Neighbor));
        for (j = 0; j < cd_args->corpus_size; j++) {
            cd_args->query_part_points[i]->neighbors[j].index = j + 1; // BIG CHANGE!!!
            cd_args->query_part_points[i]->neighbors[j].distance = distances_part_matrix[i * cd_args->corpus_size + j];
        }
    }

    for (i = 0; i < cd_args->query_part_size; i++) {
        qsort(cd_args->query_part_points[i]->neighbors, (size_t)cd_args->corpus_size, sizeof(Neighbor), qsort_compare);
        cd_args->query_part_points[i]->neighbors = realloc(cd_args->query_part_points[i]->neighbors, (size_t)cd_args->knns * sizeof(Neighbor));
    }

    // for (i = 0; i < cd_args->query_part_size; i++) {
    //     for (j = 0; j < cd_args->knns; j++) {
    //         printf("%lf ", cd_args->query_part_points[i].neighbors[j].distance);
    //     }
    //     printf("\b\n");
    // }

    free(distances_part_matrix);
    free(c_norms);
    free(q_norms);
}

void knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **set_points) {

    int i, j, t, q_part_size, split_factor;
    double **q_parts;
    pthread_t q_thread_ids[Q_SPLIT];
    calculate_distances_args **cd_args;
    Point ***set_points_parts;

    printf("k-NN search started!\n");
    
    q_part_size = (int)q_size / Q_SPLIT;
    if (q_part_size < 1) {
        q_part_size = 1;
        split_factor = q_size;
    } else {
        split_factor = Q_SPLIT;
    }

    set_points_parts = (Point ***)malloc(split_factor * sizeof(Point **));
    q_parts = (double **)malloc(split_factor * sizeof(double *));
    for (i = 0; i < split_factor; i++) {
        set_points_parts[i] = set_points + i * q_part_size;
        q_parts[i] = Q + i * d * q_part_size;
    }

    cd_args = (calculate_distances_args **)malloc(split_factor * sizeof(calculate_distances_args *));
    for (i = 0; i < split_factor; i++) {
        cd_args[i] = malloc(sizeof(calculate_distances_args));
        cd_args[i]->corpus = C;
        cd_args[i]->query_part = q_parts[i];
        cd_args[i]->corpus_size = c_size;
        cd_args[i]->query_part_size = q_part_size; // Sto teleutaio part xanetai to ypoloipo...
        cd_args[i]->query_part_points = set_points_parts[i];
        cd_args[i]->dimensions = d;
        cd_args[i]->knns = knns;
        pthread_create(&q_thread_ids[i], NULL, calculate_distances, cd_args[i]);
    }

    for (i = 0; i < split_factor; i++) {
        pthread_join(q_thread_ids[i], NULL);
    }

    printf("k-NN search finished!\n");

    // for (t = 0; t < split_factor; t++) {
    //     for (i = 0; i < cd_args[t]->query_part_size; i++) {
    //         my_idx[i + t * cd_args[t]->query_part_size] = (double *)malloc(knns * sizeof(double));
    //         my_dst[i + t * cd_args[t]->query_part_size] = (double *)malloc(knns * sizeof(double));
    //         for (j = 0; j < knns; j++) {
    //             my_idx[i + t * cd_args[t]->query_part_size][j] = (double)(cd_args[t]->query_part_points[i].neighbors[j].index + 1);
    //             my_dst[i + t * cd_args[t]->query_part_size][j] = cd_args[t]->query_part_points[i].neighbors[j].distance;
    //         }
    //     }
    // }

    free(set_points_parts);
    free(q_parts);
    for (i = 0; i < split_factor; i++) {
        free(cd_args[i]);
    }
    free(cd_args);
}
