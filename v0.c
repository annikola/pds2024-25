#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

#define MAX_SET_SPLIT 2
#define MIN_ARGS 5
#define MIN_FEATURE 1.0
#define MAX_FEATURE 100.0
#define Q_SPLIT 4
#define MAX_GB 1

typedef struct {
    double **corpus;
    double **query_part;
    double **distances;
    int corpus_size;
    int query_part_size;
    int dimensions;
} calculate_distances_args;

int qsort_compare(const void* a, const void* b);
double random_double(double min, double max);
void *calculate_distances(void *args);
double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);

int main(int argc, char *argv[]) {

    int i, j;
    int c_size, q_size, d, nns;
    int q_part_size;
    double elapsed;
    double **c, **q;
    double ***q_parts;
    double **indexes_matrix;
    struct timespec start, end;
    calculate_distances_args **cd_args;
    pthread_t q_thread_ids[Q_SPLIT];
    // const char *filename;

    /* THA FYGEI STO TELOS!!! --> */
    
    if (argc < MIN_ARGS) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    srand(time(0));

    // filename = &argv[1];
    c_size = atoi(argv[1]);
    q_size = atoi(argv[2]);
    d = atoi(argv[3]);
    nns = atoi(argv[4]);

    if (q_size * c_size * sizeof(double) / 1e9 > MAX_GB) {
        printf("Matrices are too big...look for an approximate solution!\n");
        return 0;
    }

    c = (double **)malloc(c_size * sizeof(double *));
    for (i = 0; i < c_size; i++) {
        c[i] = (double *)malloc(d * sizeof(double));
    }

    for (i = 0; i < c_size; i++) {
        for (j = 0; j < d; j++) {
            c[i][j] = random_double(MIN_FEATURE, MAX_FEATURE);
        }
    }

    q = (double **)malloc(q_size * sizeof(double *));
    for (j = 0; j < q_size; j++) {
        q[j] = (double *)malloc(d * sizeof(double));
    }

    for (i = 0; i < q_size; i++) {
        for (j = 0; j < d; j++) {
            q[i][j] = random_double(MIN_FEATURE, MAX_FEATURE);
        }
    }

    // for (i = 0; i < c_size; i++) {
    //     for (j = 0; j < d; j++) {
    //         printf("%lf ", c[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n");

    /* <-- THA FYGEI STO TELOS!!! */

    clock_gettime(CLOCK_MONOTONIC, &start);
    
    q_parts = (double ***)malloc(Q_SPLIT * sizeof(double **));
    q_part_size = (int)q_size / Q_SPLIT;
    for (i = 0; i < Q_SPLIT; i++) {
        q_parts[i] = q + i * q_part_size;
    }

    // #pragma omp parallel num_threads(Q_SPLIT)
    // {
    //     int thread_id = omp_get_thread_num();
    //     int startRow = thread_id * q_part_size;
    //     int endRow = startRow + q_part_size;
    //     calculate_distances_args *cd_args_p;

    //     cd_args = malloc(sizeof(calculate_distances_args));
    //     cd_args_p->corpus = c;
    //     cd_args_p->query_part = q;
    //     cd_args_p->corpus_size = c_size;
    //     cd_args_p->query_part_size = q_part_size; // Sto teleutaio part xanetai to ypoloipo...
    //     cd_args_p->dimensions = d;
        
    //     // Each thread processes a different chunk of rows
    //     calculate_distances(cd_args[0]);
    // }

    cd_args = (calculate_distances_args **)malloc(Q_SPLIT * sizeof(calculate_distances_args *));
    for (i = 0; i < Q_SPLIT; i++) {
        cd_args[i] = malloc(sizeof(calculate_distances_args));
        cd_args[i]->corpus = c;
        cd_args[i]->query_part = q_parts[i];
        cd_args[i]->corpus_size = c_size;
        cd_args[i]->query_part_size = q_part_size; // Sto teleutaio part xanetai to ypoloipo...
        cd_args[i]->dimensions = d;
        pthread_create(&q_thread_ids[i], NULL, calculate_distances, cd_args[i]);
    }

    for (i = 0; i < Q_SPLIT; i++) {
        pthread_join(q_thread_ids[i], NULL);
    }

    // cd_args = (calculate_distances_args **)malloc(Q_SPLIT * sizeof(calculate_distances_args *));
    // cd_args[0] = malloc(sizeof(calculate_distances_args));
    // cd_args[0]->corpus = c;
    // cd_args[0]->query_part = q;
    // cd_args[0]->corpus_size = c_size;
    // cd_args[0]->query_part_size = q_size; // Sto teleutaio part xanetai to ypoloipo...
    // cd_args[0]->dimensions = d;
    // calculate_distances(cd_args[0]);

    // for (i = 0; i < q_size; i++) {
    //     qsort(cd_args[0]->distances[i], (size_t)c_size, sizeof(double), qsort_compare);
    //     cd_args[0]->distances[i] = realloc(cd_args[0]->distances[i], nns * sizeof(double));
    // }

    clock_gettime(CLOCK_MONOTONIC, &end);

    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("kNN finished in: %lf seconds!\n\n", elapsed);

    write_2D_array_to_matfile("my_c.mat", "C", c, c_size, d);
    write_2D_array_to_matfile("my_q.mat", "Q", q, q_size, d);
    // write_2D_array_to_matfile("my_dst.mat", "ddd", cd_args[0]->distances, q_size, nns);

    free(c);
    free(q);

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

void *calculate_distances(void *args) {

    int i, j, k;
    double **distances_part_matrix;
    calculate_distances_args *cd_args;

    cd_args = (calculate_distances_args *)args;

    distances_part_matrix = (double **)malloc((cd_args->query_part_size) * sizeof(double *));
    for (k = 0; k < cd_args->query_part_size; k++) {
        distances_part_matrix[k] = (double *)malloc(cd_args->corpus_size * sizeof(double));
        for (i = 0; i < cd_args->corpus_size; i++) {
            distances_part_matrix[k][i] = 0.0;
            for (j = 0; j < cd_args->dimensions; j++) {
                distances_part_matrix[k][i] += pow(cd_args->corpus[i][j] - cd_args->query_part[k][j], 2);
            }
            distances_part_matrix[k][i] = sqrt(distances_part_matrix[k][i]);
        }
    }

    cd_args->distances = distances_part_matrix;

    // for (i = 0; i < cd_args->query_part_size; i++) {
    //     for (j = 0; j < cd_args->corpus_size; j++) {
    //         printf("%lf ", distances_part_matrix[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n");
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
