#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <cblas.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

#define MAX_SET_SPLIT 2
#define MIN_ARGS 5
#define MIN_FEATURE 1.0
#define MAX_FEATURE 100.0
#define Q_SPLIT 4
#define MAX_GB 1

typedef struct {
    double *corpus;
    double *query_part;
    double *distances;
    int corpus_size;
    int query_part_size;
    int dimensions;
    int knns;
} calculate_distances_args;

int qsort_compare(const void* a, const void* b);
double random_double(double min, double max);
double *calculate_norms(double *matrix2D, int m_size, int d);
void *calculate_distances(void *args);
double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);

int main(int argc, char *argv[]) {

    int i, j;
    int c_size, q_size, d, q_part_size, knns;
    double elapsed;
    double **q_parts;
    double **indexes_matrix;
    struct timespec start, end;
    double *C, *Q, *result;
    calculate_distances_args **cd_args;
    pthread_t q_thread_ids[Q_SPLIT];
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

    Q = (double *)malloc(q_size * d * sizeof(double));
    for (i = 0; i < q_size; i++) {
        for (j = 0; j < d; j++) {
            Q[i * d + j] = random_double(MIN_FEATURE, MAX_FEATURE);
        }
    }
    
    result = (double *)malloc(c_size * q_size * sizeof(double));

    /* <-- THA FYGEI STO TELOS!!! */

    clock_gettime(CLOCK_MONOTONIC, &start);
    
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
        cd_args[i]->dimensions = d;
        cd_args[i]->knns = knns;
        pthread_create(&q_thread_ids[i], NULL, calculate_distances, cd_args[i]);
    }

    for (i = 0; i < Q_SPLIT; i++) {
        pthread_join(q_thread_ids[i], NULL);
    }

    // clock_gettime(CLOCK_MONOTONIC, &start);

    // cd_args = (calculate_distances_args **)malloc(Q_SPLIT * sizeof(calculate_distances_args *));
    // cd_args[0] = malloc(sizeof(calculate_distances_args));
    // cd_args[0]->corpus = c;
    // cd_args[0]->query_part = q;
    // cd_args[0]->corpus_size = c_size;
    // cd_args[0]->query_part_size = q_size; // Sto teleutaio part xanetai to ypoloipo...
    // cd_args[0]->dimensions = d;
    // calculate_distances(cd_args[0]);

    // for (i = 0; i < Q_SPLIT; i++) {
    //     for (j = 0; j < cd_args[i]->query_part_size; j++) {
    //         qsort(cd_args[i]->distances[j], (size_t)c_size, sizeof(double), qsort_compare);
    //         cd_args[i]->distances[j] = realloc(cd_args[i]->distances[j], nns * sizeof(double));
    //     }
    // }

    // clock_gettime(CLOCK_MONOTONIC, &end);

    // elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    // printf("kNN finished in: %lf seconds!\n\n", elapsed);

    // clock_gettime(CLOCK_MONOTONIC, &start);

    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, c_size, c_size, d, 1.0, C, d, C, d, 0.0, result, c_size);

    clock_gettime(CLOCK_MONOTONIC, &end);

    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("openBLAS finished in: %lf seconds!\n\n", elapsed);

    // for (int t = 0; t < Q_SPLIT; t++) {
    //     printf("Result matrix result:\n");
    //     for (i = 0; i < cd_args[t]->corpus_size; i++) {
    //         for (j = 0; j < cd_args[t]->query_part_size; j++) {
    //             printf("%lf ", cd_args[t]->distances[i * cd_args[t]->query_part_size + j]);
    //         }
    //         printf("\n");
    //     }
    // }

    // printf("Result matrix result:\n");
    // for (i = 0; i < c_size; i++) {
    //     for (j = 0; j < c_size; j++) {
    //         result[i * c_size + j] = random_double(MIN_FEATURE, MAX_FEATURE);
    //         printf("%f ", result[i * c_size + j]);
    //     }
    //     printf("\n");
    // }

    // write_2D_array_to_matfile("my_c.mat", "C", c, c_size, d);
    // write_2D_array_to_matfile("my_q.mat", "Q", q, q_size, d);
    // write_2D_array_to_matfile("my_dst.mat", "ddd", cd_args[0]->distances, q_size, nns);

    // free(C);
    // free(Q);

    return 0;
}

int qsort_compare(const void *a, const void *b) {

    double diff = *(double *)a - *(double *)b;

    if (diff < 0) {
        return -1;
    } else if (diff > 0) {
        return 1;
    } else {
        return 0;
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
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, cd_args->corpus_size, cd_args->query_part_size, cd_args->dimensions, -2.0, cd_args->corpus, cd_args->dimensions, cd_args->query_part, cd_args->dimensions, 0.0, matrixmul, cd_args->query_part_size);
    
    c_norms = calculate_norms(cd_args->corpus, cd_args->corpus_size, cd_args->dimensions);
    q_norms = calculate_norms(cd_args->query_part, cd_args->query_part_size, cd_args->dimensions);

    distances_part_matrix = (double *)malloc(cd_args->corpus_size * cd_args->query_part_size * sizeof(double));
    for (i = 0; i < cd_args->corpus_size; ++i) {
        for (j = 0; j < cd_args->query_part_size; ++j) {
            distances_part_matrix[i * cd_args->query_part_size + j] = c_norms[i] + q_norms[j] + matrixmul[i * cd_args->query_part_size + j];
            // Take the square root to get the Euclidean distance
            distances_part_matrix[i * cd_args->query_part_size + j] = sqrt(distances_part_matrix[i * cd_args->query_part_size + j]);
        }
    }

    // qsort(distances_part_matrix, (size_t)cd_args->corpus_size, sizeof(double), qsort_compare);
    // cd_args->distances = realloc(distances_part_matrix, cd_args->knns * sizeof(double));

    cd_args->distances = distances_part_matrix;
}

double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d) {

    MATFile *pmat;
    mxArray *array_ptr;
    size_t i, j;
    double *data;
    double **c;
    const char *varname = "C"; // Name of the variable to read

    pmat = matOpen(filename, "r");
    if (pmat == NULL) {
        fprintf(stderr, "Error opening file data.mat\n");
        return NULL;
    }

    // Read a variable from the .mat file
    array_ptr = matGetVariable(pmat, varname);
    if (array_ptr == NULL) {
        fprintf(stderr, "Error reading variable %s from file\n", varname);
        matClose(pmat);
        return NULL;
    }

    // Check if the variable is of the expected type (for example, double matrix)
    if (mxIsDouble(array_ptr) && !mxIsComplex(array_ptr)) {
        // Get the data pointer and array dimensions
        data = mxGetPr(array_ptr);
        size_t rows = mxGetM(array_ptr);
        size_t cols = mxGetN(array_ptr);

        *c_size = rows;
        *d = cols;

        c = (double **)malloc(rows * sizeof(double *));
        for (i = 0; i < rows; i++) {
            c[i] = (double *)malloc(cols * sizeof(double));
        }

        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                c[i][j] = data[i + j * rows];
            }
        }
    } else {
        fprintf(stderr, "Variable %s is not a double matrix\n", varname);
        return NULL;
    }

    // Clean up
    mxDestroyArray(array_ptr); // Free the mxArray
    matClose(pmat);            // Close the MAT-file

    return c;
}

void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d) {

    int i, j;
    double *matData;
    MATFile *matFile;
    mxArray *matArray;

    matFile = matOpen(filename, "w");
    if (matFile == NULL) {
        printf("Error opening MAT-file %s\n", filename);
        return;
    }

    matArray = mxCreateDoubleMatrix(c_size, d, mxREAL);
    if (matArray == NULL) {
        printf("Could not create mxArray.\n");
        matClose(matFile);
        return;
    }

    // matData = mxGetPr(matArray);
    // for (i = 0; i < c_size; i++) {
    //     memcpy(matData + i * d, _2D_array[i], d * sizeof(double));
    // }

    // MATLAB saves in column major order, therefore we take the transpose...
    matData = mxGetPr(matArray);
    for (i = 0; i < c_size; i++) {
        for (j = 0; j < d; j++) {
            matData[j * c_size + i] = _2D_array[i][j];
        }
    }

    matPutVariable(matFile, array_name, matArray);

    mxDestroyArray(matArray);
    matClose(matFile);
    printf("2D array written to %s successfully.\n", filename);
}
