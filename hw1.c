#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

#define MAX_THREADS 5
#define MAX_SET_SPLIT 2
#define MIN_ARGS 5

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

typedef struct {
    struct hyper_set *root;
    struct hyper_set *leaf;
    double *q_point;
} hyper_traverse_args;

void *hyper_binary_traverse(void *traverse_args);
struct hyper_set *hyper_search(struct hyper_set *hyper_subset, double *q_point);
void *hyper_binary_split(void *hyper_subset);
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
    struct hyper_set **root_hyper_sets;
    hyper_traverse_args **traverse_args;

    /* THA FYGEI STO TELOS!!! --> */
    
    if (argc < MIN_ARGS) {
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

    start_t = clock();
    root_hyper_sets = (struct hyper_set **)malloc(MAX_THREADS * sizeof(struct hyper_set *));
    for (t = 0; t < MAX_THREADS; t++) {
        root_hyper_sets[t] = malloc(sizeof(struct hyper_set));
        root_hyper_sets[t]->indexes = NULL;
        root_hyper_sets[t]->neighbors = c;
        root_hyper_sets[t]->set_size = c_size;
        root_hyper_sets[t]->d = d;
        root_hyper_sets[t]->coeffs = NULL; // Technically a zero-vector...
        root_hyper_sets[t]->b_value = 0.0;
        root_hyper_sets[t]->depth = depth;
        root_hyper_sets[t]->f_ptr = NULL;
        pthread_create(&thread_ids[t], NULL, hyper_binary_split, root_hyper_sets[t]);
    }

    for (t = 0; t < MAX_THREADS; t++) {
        pthread_join(thread_ids[t], NULL);
    }

    printf("Multithreaded application finished in: %lf seconds!\n\n", (double)(clock() - start_t) / CLOCKS_PER_SEC);

    printf("Point of interest: ");
    for (j = 0; j < d; j++) {
        printf("%lf ", q[0][j]);
    }
    printf("\b\n\n");
    
    traverse_args = (hyper_traverse_args **)malloc(MAX_THREADS * sizeof(hyper_traverse_args *));
    for (t = 0; t < MAX_THREADS; t++) {
        traverse_args[t] = malloc(sizeof(hyper_binary_traverse));
        traverse_args[t]->root = (struct hyper_set *)root_hyper_sets[t];
        traverse_args[t]->leaf = NULL;
        traverse_args[t]->q_point = q[0];
        pthread_create(&thread_ids[t], NULL, hyper_binary_traverse, traverse_args[t]);
    }

    for (t = 0; t < MAX_THREADS; t++) {
        pthread_join(thread_ids[t], NULL);
    }
    printf("Tree scanning finished in: %lf seconds!\n\n", (double)(clock() - start_t) / CLOCKS_PER_SEC);

    // Print all leafs!
    for (t = 0; t < MAX_THREADS; t++) {
        for (i = 0; i < traverse_args[t]->leaf->set_size; i++) {
            for (j = 0; j < traverse_args[t]->leaf->d; j++) {
                printf("%lf ", traverse_args[t]->leaf->neighbors[i][j]);
            }
            printf("\b\n");
        }
        printf("\n\n");
    }

    free(c);
    free(q);
    free(root_hyper_sets);
    free(traverse_args);

    return 0;
}

void *hyper_binary_traverse(void *traverse_args) {

    hyper_traverse_args *curr_args;

    curr_args = (hyper_traverse_args *) traverse_args;

    curr_args->leaf = hyper_search(curr_args->root, curr_args->q_point);
    while ((curr_args->leaf->f_ptr) != NULL) {
        curr_args->leaf = hyper_search(curr_args->leaf, curr_args->q_point);
    }

}

struct hyper_set *hyper_search(struct hyper_set *hyper_subset, double *q_point) {

    int i, j;
    double hyper_position;

    if (hyper_subset->coeffs == NULL) {
        printf("YES\n");
    }

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

void *hyper_binary_split(void *hyper_subset_void) {

    int i, j;
    double *random_point_1, *random_point_2;
    double *midpoint, *normal_vector;
    double beta, hyper_position;
    struct hyper_set *hyper_subset, *new_hyper_subset_1, *new_hyper_subset_2;

    hyper_subset = (struct hyper_set *)hyper_subset_void;

    random_point_1 = hyper_subset->neighbors[rand() % hyper_subset->set_size];
    random_point_2 = hyper_subset->neighbors[rand() % hyper_subset->set_size];

    // RARE CASE!!!
    while (is_duplicate(random_point_1, random_point_2, hyper_subset->d)) {
        random_point_2 = hyper_subset->neighbors[rand() % hyper_subset->set_size];
    }

    midpoint = (double *)malloc(hyper_subset->d * sizeof(double));
    for (j = 0; j < hyper_subset->d; j++) {
        midpoint[j] = (random_point_1[j] + random_point_2[j]) / 2.0;
    }

    normal_vector = (double *)malloc(hyper_subset->d * sizeof(double));
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

    printf("Thread finished!\n");
}
