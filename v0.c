#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <cblas.h>

#define MAX_SET_SPLIT 2
#define MIN_ARGS 5
#define MIN_FEATURE 1.0
#define MAX_FEATURE 100.0
#define Q_SPLIT 4
#define MAX_GB 1

typedef struct {
    double distance;
    int index;
} Neighbor;

typedef struct {
    Neighbor *neighbors;
    int index;
} Point;

typedef struct {
    double *corpus;
    double *query_part;
    Point *query_part_points;
    int corpus_size;
    int query_part_size;
    int dimensions;
    int knns;
} calculate_distances_args;

int qsort_compare(const void* a, const void* b);
double random_double(double min, double max);
double *calculate_norms(double *matrix2D, int m_size, int d);
void *calculate_distances(void *args);
calculate_distances_args **knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns);
double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);

int main(int argc, char *argv[]) {

    int i, j, t;
    int c_size, q_size, d, knns;
    // int **my_idx;
    double elapsed;
    double *C, *Q;
    double **my_c, **my_q, **my_idx, **my_dst;
    struct timespec start, end;
    calculate_distances_args **cd_args;
    // const char *filename;
    
    if (argc < MIN_ARGS) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    /* THA FYGEI STO TELOS!!! --> */

    srand(time(0));

    // filename = &argv[1];
    c_size = atoi(argv[1]);
    q_size = atoi(argv[2]);
    d = atoi(argv[3]);
    knns = atoi(argv[4]);

    if (q_size * c_size * sizeof(double) / 10e9 > MAX_GB) {
        printf("Matrices are too big...look for an approximate solution!\n");
        return 0;
    }

    C = (double *)malloc(c_size * d * sizeof(double));
    for (i = 0; i < c_size; i++) {
        for (j = 0; j < d; j++) {
            C[i * d + j] = random_double(MIN_FEATURE, MAX_FEATURE);
        }
    }

    // This is the format that the matlab_write accepts...
    my_c = (double **)malloc(c_size * sizeof(double *));
    for (i = 0; i < c_size; i++) {
        my_c[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            my_c[i][j] = C[i * d + j];
        }
    }

    Q = (double *)malloc(q_size * d * sizeof(double));
    for (i = 0; i < q_size; i++) {
        for (j = 0; j < d; j++) {
            Q[i * d + j] = random_double(MIN_FEATURE, MAX_FEATURE);
        }
    }

    my_q = (double **)malloc(q_size * sizeof(double *));
    for (i = 0; i < q_size; i++) {
        my_q[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            my_q[i][j] = Q[i * d + j];
        }
    }

    /* <-- THA FYGEI STO TELOS!!! */

    clock_gettime(CLOCK_MONOTONIC, &start);

    cd_args = knn_search(C, Q, c_size, q_size, d, knns);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("k-NN search finished in: %lf seconds!\n\n", elapsed);

    // my_idx = (int **)malloc(q_size * sizeof(int *));
    // for (int t = 0; t < Q_SPLIT; t++) {
    //     for (i = 0; i < q_part_size; i++) {
    //         my_idx[i + t * q_part_size] = (int *)malloc(knns * sizeof(int));
    //         for (j = 0; j < knns; j++) {
    //             my_idx[i + t * q_part_size][j] = cd_args[t]->query_part_points[i].neighbors[j].index;
    //         }
    //     }
    // }

    my_idx = (double **)malloc(q_size * sizeof(double *));
    my_dst = (double **)malloc(q_size * sizeof(double *));
    for (t = 0; t < Q_SPLIT; t++) {
        for (i = 0; i < cd_args[t]->query_part_size; i++) {
            my_idx[i + t * cd_args[t]->query_part_size] = (double *)malloc(knns * sizeof(double));
            my_dst[i + t * cd_args[t]->query_part_size] = (double *)malloc(knns * sizeof(double));
            for (j = 0; j < knns; j++) {
                my_idx[i + t * cd_args[t]->query_part_size][j] = (double)(cd_args[t]->query_part_points[i].neighbors[j].index + 1);
                my_dst[i + t * cd_args[t]->query_part_size][j] = cd_args[t]->query_part_points[i].neighbors[j].distance;
            }
        }
    }

    printf("Writing to mat files...\n");

    write_2D_array_to_matfile("my_c.mat", "C", my_c, c_size, d);
    write_2D_array_to_matfile("my_q.mat", "Q", my_q, q_size, d);
    write_2D_array_to_matfile("my_idx.mat", "iii", my_idx, q_size, knns);
    write_2D_array_to_matfile("my_dst.mat", "ddd", my_dst, q_size, knns);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Total time: %lf seconds!\n", elapsed);

    free(C);
    free(Q);
    free(my_idx);
    free(my_dst);
    free(my_c);
    free(my_q);

    return 0;
}

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

double random_double(double min, double max) {

    double scale;

    scale = rand() / (double) RAND_MAX;

    return min + scale * (max - min);
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
        cd_args->query_part_points[i].neighbors = (Neighbor *)malloc(cd_args->corpus_size * sizeof(Neighbor));
        for (j = 0; j < cd_args->corpus_size; j++) {
            cd_args->query_part_points[i].neighbors[j].index = j;
            cd_args->query_part_points[i].neighbors[j].distance = distances_part_matrix[i * cd_args->corpus_size + j];
        }
    }

    for (i = 0; i < cd_args->query_part_size; i++) {
        qsort(cd_args->query_part_points[i].neighbors, (size_t)cd_args->corpus_size, sizeof(Neighbor), qsort_compare);
        cd_args->query_part_points[i].neighbors = realloc(cd_args->query_part_points[i].neighbors, (size_t)cd_args->knns * sizeof(Neighbor));
    }

    // for (i = 0; i < cd_args->query_part_size; i++) {
    //     for (j = 0; j < cd_args->knns; j++) {
    //         printf("%lf ", cd_args->query_part_points[i].neighbors[j].distance);
    //     }
    //     printf("\b\n");
    // }
}

calculate_distances_args **knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns) {

    int i, q_part_size;
    double **q_parts;
    pthread_t q_thread_ids[Q_SPLIT];
    calculate_distances_args **cd_args;

    printf("k-NN search started!\n");
    
    q_parts = (double **)malloc(Q_SPLIT * sizeof(double *));
    q_part_size = (int)q_size / Q_SPLIT;
    for (i = 0; i < Q_SPLIT; i++) {
        q_parts[i] = Q + i * d * q_part_size;
    }

    cd_args = (calculate_distances_args **)malloc(Q_SPLIT * sizeof(calculate_distances_args *));
    for (i = 0; i < Q_SPLIT; i++) {
        cd_args[i] = malloc(sizeof(calculate_distances_args));
        cd_args[i]->corpus = C;
        cd_args[i]->query_part = q_parts[i];
        cd_args[i]->corpus_size = c_size;
        cd_args[i]->query_part_size = q_part_size; // Sto teleutaio part xanetai to ypoloipo...
        cd_args[i]->query_part_points = (Point *)malloc(q_part_size * sizeof(Point));
        cd_args[i]->dimensions = d;
        cd_args[i]->knns = knns;
        pthread_create(&q_thread_ids[i], NULL, calculate_distances, cd_args[i]);
    }

    for (i = 0; i < Q_SPLIT; i++) {
        pthread_join(q_thread_ids[i], NULL);
    }

    return cd_args;
}
