#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

#define MAX_THREADS 5
#define MAX_SET_SPLIT 2

typedef struct {
    double *q_curr;
    int ccp_length;
    int d;
} ThreadArgs;

struct hyper_set {
    int *indexes;
    double **neighbors;
    int set_size;
    int d;
    double *coeffs;
    double b_value;
    int depth;
    struct hyper_set **f_ptr;
};

struct hyper_set *hyper_search(struct hyper_set *hyper_subset, double *q_point);
void *hyper_binary_split(struct hyper_set *hyper_subset);
int is_duplicate(double *point1, double *point2, int d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);
double random_double(double min, double max);
void *distance_calculator(void *args);

int main(int argc, char *argv[]) {

    int i, j, t;
    int c_size, q_size, d, depth;
    double **c, **q;
    pthread_t thread_ids[MAX_THREADS];
    clock_t start_t;
    ThreadArgs *args;
    struct hyper_set *root_hyper_set, *temp;

    /* THA FYGEI STO TELOS!!! --> */
    
    if (argc < 5) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    srand(time(0));

    c_size = atoi(argv[1]);
    q_size = atoi(argv[2]);
    d = atoi(argv[3]);
    depth = atoi(argv[4]);

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

    write_2D_array_to_matfile("test.mat", "C", c, c_size, d);
    write_2D_array_to_matfile("test2.mat", "Q", q, q_size, d);
    printf("\n");

    for (i = 0; i < c_size; i++) {
        for (j = 0; j < d; j++) {
            printf("%lf ", c[i][j]);
        }
        printf("\b\n");
    }
    printf("\n");

    /* <-- THA FYGEI STO TELOS!!! */

    root_hyper_set = (struct hyper_set *)malloc(sizeof(struct hyper_set));
    root_hyper_set->indexes = NULL;
    root_hyper_set->neighbors = c;
    root_hyper_set->set_size = c_size;
    root_hyper_set->d = d;
    root_hyper_set->coeffs = NULL; // Technically a zero-vector...
    root_hyper_set->b_value = 0.0;
    root_hyper_set->depth = depth;
    root_hyper_set->f_ptr = NULL;

    start_t = clock();
    hyper_binary_split(root_hyper_set);
    printf("Multithreaded application finished in: %lf seconds!\n\n", (double)(clock() - start_t) / CLOCKS_PER_SEC);

    for (j = 0; j < d; j++) {
        printf("%lf ", q[0][j]);
    }
    printf("\b\n\n");
    
    temp = hyper_search(root_hyper_set, q[0]);
    while ((temp->f_ptr) != NULL) {
        temp = hyper_search(temp, q[0]);
    }

    // temp points to a leaf!!! (The leaf is ONLY traversed below...)
    for (i = 0; i < temp->set_size; i++) {
        for (j = 0; j < temp->d; j++) {
            printf("%lf ", temp->neighbors[i][j]);
        }
        printf("\b\n");
    }
    printf("\n\n");

    printf("Tree scanning finished in: %lf seconds!\n", (double)(clock() - start_t) / CLOCKS_PER_SEC);

    // for (k = 0; k < MAX_CLUSTERS; k++) {
    //     for (i = 0; i < temp[k] - clusters[k].neighbors; i++) {
    //         for (j = 0; j < d; j++) {
    //             printf("%lf ", clusters[k].neighbors[i][j]);
    //         }
    //         printf("\b\n");
    //     }
    //     printf("\n");
    // }

    // start_t = clock();
    // args = (ThreadArgs *)malloc(sizeof(ThreadArgs)); ???
    // args = malloc(sizeof(ThreadArgs));
    // for (t = 0; t < MAX_THREADS; t++) {
    //     args->q_curr = q[0];
    //     args->ccp_length = 0 ???
    //     args->d = d;
    //     pthread_create(&thread_ids[t], NULL, distance_calculator, args);
    // }

    // for (t = 0; t < MAX_THREADS; t++) {
    //     pthread_join(thread_ids[t], NULL);
    // }

    // printf("Multithreaded application finished in: %lf seconds!\n", (double)(clock() - start_t) / CLOCKS_PER_SEC);

    free(c);
    free(q);
    free(root_hyper_set);

    return 0;
}

struct hyper_set *hyper_search(struct hyper_set *hyper_subset, double *q_point) {

    int i, j;
    double hyper_position;

    hyper_position = 0.0;
    for (j = 0; j < hyper_subset->d; j++) {
        hyper_position += hyper_subset->coeffs[j] * q_point[j];
    }

    // for (i = 0; i < hyper_subset->f_ptr[0]->set_size; i++) {
    //     for (j = 0; j < hyper_subset->d; j++) {
    //         printf("%lf ", hyper_subset->f_ptr[0]->neighbors[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n\n");

    // for (i = 0; i < hyper_subset->f_ptr[1]->set_size; i++) {
    //     for (j = 0; j < hyper_subset->d; j++) {
    //         printf("%lf ", hyper_subset->f_ptr[1]->neighbors[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n\n");

    if (hyper_position > hyper_subset->b_value) {
        return hyper_subset->f_ptr[0];
    } else {
        return hyper_subset->f_ptr[1];
    }
}

void *hyper_binary_split(struct hyper_set *hyper_subset) {

    int i, j;
    double *random_point_1, *random_point_2;
    double midpoint[hyper_subset->d], normal_vector[hyper_subset->d];
    double beta, hyper_position;
    struct hyper_set *new_hyper_subset_1, *new_hyper_subset_2;

    random_point_1 = hyper_subset->neighbors[rand() % hyper_subset->set_size];
    random_point_2 = hyper_subset->neighbors[rand() % hyper_subset->set_size];

    // RARE CASE!!!
    while (is_duplicate(random_point_1, random_point_2, hyper_subset->d)) {
        random_point_2 = hyper_subset->neighbors[rand() % hyper_subset->set_size];
    }

    for (j = 0; j < hyper_subset->d; j++) {
        midpoint[j] = (random_point_1[j] + random_point_2[j]) / 2.0;
    }

    for (j = 0; j < hyper_subset->d; j++) {
        normal_vector[j] = random_point_2[j] - random_point_1[j];
    }

    beta = 0.0;
    for (j = 0; j < hyper_subset->d; j++) {
        beta += normal_vector[j] * midpoint[j];
    }

    hyper_subset->coeffs = normal_vector;
    hyper_subset->b_value = beta;

    new_hyper_subset_1 = (struct hyper_set *)malloc(sizeof(struct hyper_set));
    new_hyper_subset_1->indexes = (int *)malloc(hyper_subset->set_size * sizeof(int));
    new_hyper_subset_1->neighbors = (double **)malloc(hyper_subset->set_size * sizeof(double *));

    new_hyper_subset_2 = (struct hyper_set *)malloc(sizeof(struct hyper_set));
    new_hyper_subset_2->indexes = (int *)malloc(hyper_subset->set_size * sizeof(int));
    new_hyper_subset_2->neighbors = (double **)malloc(hyper_subset->set_size * sizeof(double *));

    new_hyper_subset_1->depth = new_hyper_subset_2->depth = hyper_subset->depth;
    new_hyper_subset_1->d = new_hyper_subset_2->d = hyper_subset->d;
    new_hyper_subset_1->set_size = new_hyper_subset_2->set_size = 0;
    for (i = 0; i < hyper_subset->set_size; i++) {
        hyper_position = 0.0;
        for (j = 0; j < hyper_subset->d; j++) {
            hyper_position += normal_vector[j] * hyper_subset->neighbors[i][j];
        }

        if (hyper_position > beta) {
            new_hyper_subset_1->indexes[new_hyper_subset_1->set_size] = i;
            new_hyper_subset_1->neighbors[new_hyper_subset_1->set_size] = hyper_subset->neighbors[i];
            new_hyper_subset_1->set_size++;
        }
        else {
            new_hyper_subset_2->indexes[new_hyper_subset_2->set_size] = i;
            new_hyper_subset_2->neighbors[new_hyper_subset_2->set_size] = hyper_subset->neighbors[i];
            new_hyper_subset_2->set_size++;
        }
    }

    new_hyper_subset_1->indexes = realloc(new_hyper_subset_1->indexes, hyper_subset->set_size * sizeof(int));
    new_hyper_subset_1->neighbors = realloc(new_hyper_subset_1->neighbors, new_hyper_subset_1->set_size * sizeof(double *));

    new_hyper_subset_2->indexes = realloc(new_hyper_subset_2->indexes, hyper_subset->set_size * sizeof(int));
    new_hyper_subset_2->neighbors = realloc(new_hyper_subset_2->neighbors, new_hyper_subset_2->set_size * sizeof(double *));

    hyper_subset->f_ptr = (struct hyper_set **)malloc(2 * sizeof(struct hyper_set *));
    hyper_subset->f_ptr[0] = new_hyper_subset_1;
    hyper_subset->f_ptr[1] = new_hyper_subset_2;

    if (new_hyper_subset_1->set_size > hyper_subset->depth) {
        hyper_binary_split(new_hyper_subset_1);
    } else {
        new_hyper_subset_1->f_ptr = NULL;
    }
    
    if (new_hyper_subset_2->set_size > hyper_subset->depth) {
        hyper_binary_split(new_hyper_subset_2);
    } else {
        new_hyper_subset_2->f_ptr = NULL;
    }

}

int is_duplicate(double *point1, double *point2, int d) {

    int j;

    for (j = 0; j < d; j++) {
        if (point1[j] != point2[j]) {
            return 0;
        }
    }

    return 1;
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

double random_double(double min, double max) {

    double scale;

    scale = rand() / (double) RAND_MAX;

    return min + scale * (max - min);
}

void *distance_calculator(void *args) {

    int i, j;
    ThreadArgs *cargs;

    cargs = (ThreadArgs *)args;

    printf("Thread finished!\n");
}
