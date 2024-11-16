#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <cblas.h>
#include <unistd.h>

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
    Point **corpus_points;
    Point **query_part_points;
    int corpus_size;
    int query_part_size;
    int dimensions;
    int knns;
    int stitching;
} calculate_distances_args;

int qsort_compare(const void* a, const void* b);
int is_a_neighbor(Neighbor *neighbors, int knns, int current_index);
double *calculate_norms(double *matrix2D, int m_size, int d);
void *calculate_distances(void *args);
void knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **set_points, int stitching);

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

int is_a_neighbor(Neighbor *neighbors, int knns, int current_index) {
    
    int i;

    for (i = 0; i < knns; i++) {
        if (neighbors[i].index == current_index) {
            return 1;
        }
    }

    return 0;
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

    int i, j;
    double curr_distance;
    double *c_norms, *q_norms, *matrixmul;
    calculate_distances_args *cd_args;
    Neighbor *curr_neighbors;

    cd_args = (calculate_distances_args *)args;

    matrixmul = (double *)malloc(cd_args->corpus_size * cd_args->query_part_size * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, cd_args->query_part_size, cd_args->corpus_size, cd_args->dimensions, -2.0, cd_args->query_part, cd_args->dimensions, cd_args->corpus, cd_args->dimensions, 0.0, matrixmul, cd_args->corpus_size);

    c_norms = calculate_norms(cd_args->corpus, cd_args->corpus_size, cd_args->dimensions);
    q_norms = calculate_norms(cd_args->query_part, cd_args->query_part_size, cd_args->dimensions);

    if (!cd_args->stitching) {
        for (i = 0; i < cd_args->query_part_size; i++) {
            cd_args->query_part_points[i]->neighbors = (Neighbor *)malloc(cd_args->corpus_size * sizeof(Neighbor));
            for (j = 0; j < cd_args->corpus_size; j++) {
                cd_args->query_part_points[i]->neighbors[j].index = cd_args->corpus_points[j]->index; // BIG CHANGE!!!

                curr_distance = q_norms[i] + c_norms[j] + matrixmul[i * cd_args->corpus_size + j];
                if (curr_distance < 0) {
                    curr_distance = 0.0;
                } else {
                    curr_distance = sqrt(curr_distance);
                }
                cd_args->query_part_points[i]->neighbors[j].distance = curr_distance;
            }
        }

        for (i = 0; i < cd_args->query_part_size; i++) {
            qsort(cd_args->query_part_points[i]->neighbors, (size_t)cd_args->corpus_size, sizeof(Neighbor), qsort_compare);
            curr_neighbors = cd_args->query_part_points[i]->neighbors;
            cd_args->query_part_points[i]->neighbors = (Neighbor *)malloc(cd_args->knns * sizeof(Neighbor));
            memcpy(cd_args->query_part_points[i]->neighbors, curr_neighbors, cd_args->knns * sizeof(Neighbor));
            // cd_args->query_part_points[i]->neighbors = realloc(cd_args->query_part_points[i]->neighbors, cd_args->knns * sizeof(Neighbor));
            free(curr_neighbors);
        }
    } else {
        for (i = 0; i < cd_args->query_part_size; i++) {
            for (j = 0; j < cd_args->corpus_size; j++) {
                curr_distance = q_norms[i] + c_norms[j] + matrixmul[i * cd_args->corpus_size + j];
                if (curr_distance < 0) {
                    curr_distance = 0.0;
                } else {
                    curr_distance = sqrt(curr_distance);
                }

                if (!is_a_neighbor(cd_args->query_part_points[i]->neighbors, cd_args->knns, cd_args->corpus_points[j]->index) && (curr_distance < cd_args->query_part_points[i]->neighbors[cd_args->knns - 1].distance)) {
                    cd_args->query_part_points[i]->neighbors[cd_args->knns - 1].index = cd_args->corpus_points[j]->index;
                    cd_args->query_part_points[i]->neighbors[cd_args->knns - 1].distance = curr_distance;
                    qsort(cd_args->query_part_points[i]->neighbors, (size_t)cd_args->knns, sizeof(Neighbor), qsort_compare);
                }
            }
        }
    }

    free(matrixmul);
    free(c_norms);
    free(q_norms);
}

void knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **set_points, int stitching) {

    int i, j, q_part_size, split_factor;
    double **q_parts;
    pthread_t q_thread_ids[Q_SPLIT];
    calculate_distances_args **cd_args;
    Point ***set_points_parts;

    // THEORITIKA DEN EINAI TELEIO GIA C != Q
    if (c_size <= 1 && q_size <= 1) {
        return;
    }

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
        if (i == split_factor - 1 && q_size > Q_SPLIT) {
            cd_args[i]->query_part_size = q_part_size + q_size % Q_SPLIT; // Sto teleutaio part prostithetai kai to ypoloipo...
        } else {
            cd_args[i]->query_part_size = q_part_size;
        }
        cd_args[i]->corpus_points = set_points;
        cd_args[i]->query_part_points = set_points_parts[i];
        cd_args[i]->dimensions = d;
        cd_args[i]->knns = knns;
        cd_args[i]->stitching = stitching;
        pthread_create(&q_thread_ids[i], NULL, calculate_distances, cd_args[i]);
    }

    for (i = 0; i < split_factor; i++) {
        pthread_join(q_thread_ids[i], NULL);
    }

    printf("k-NN search finished!\n");

    free(set_points_parts);
    free(q_parts);
    for (i = 0; i < split_factor; i++) {
        free(cd_args[i]);
    }
    free(cd_args);
}
