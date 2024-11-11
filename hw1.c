#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

#define MAX_THREADS 4
#define MAX_SET_SPLIT 2
#define MIN_ARGS 2

struct hyper_set {
    double **neighbors;
    struct hyper_set **f_ptr;
    double *coeffs;
    int *indexes;
    double b_value;
    int set_size;
    int d;
    int depth;
};

typedef struct {
    struct hyper_set *root;
    struct hyper_set *leaf;
    double *q_point;
} hyper_traverse_args;

int compare_rows(const void *a, const void *b);
void *hyper_binary_traverse(void *traverse_args);
struct hyper_set *hyper_search(struct hyper_set *hyper_subset, double *q_point);
void *hyper_binary_split(void *hyper_subset);
int is_duplicate(double *point1, double *point2, int d);
void *distance_calculator(void *args);
double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);

int main(int argc, char *argv[]) {

    int i, j, t, leafs_size, depth;
    size_t c_size, d;
    int **leafs;
    // size_t c_size, d;
    double elapsed;
    double **c, **q;
    pthread_t thread_ids[MAX_THREADS];
    clock_t start_t;
    struct timespec start, end;
    struct hyper_set **root_hyper_sets;
    hyper_traverse_args **traverse_args;
    const char *filename;

    /* THA FYGEI STO TELOS!!! --> */
    
    if (argc < MIN_ARGS) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    filename = argv[1];
    depth = atoi(argv[2]);

    c = read_2D_array_from_matfile("big_set.mat", &c_size, &d);

    // for (i = 0; i < c_size; i++) {
    //     for (j = 0; j < d; j++) {
    //         printf("%lf ", c[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n");

    /* <-- THA FYGEI STO TELOS!!! */

    clock_gettime(CLOCK_MONOTONIC, &start);
    // start_t = clock();
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

    clock_gettime(CLOCK_MONOTONIC, &end);

    // Calculate the elapsed time in seconds
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("Multithreaded application finished in: %lf seconds!\n\n", elapsed);
    // printf("Multithreaded application finished in: %lf seconds!\n\n", (double)(clock() - start_t) / CLOCKS_PER_SEC);
    
    // for (i = 0; i < 10000; i++) {
    //     traverse_args = (hyper_traverse_args **)malloc(MAX_THREADS * sizeof(hyper_traverse_args *));
    //     for (t = 0; t < MAX_THREADS; t++) {
    //         traverse_args[t] = malloc(sizeof(hyper_binary_traverse));
    //         traverse_args[t]->root = (struct hyper_set *)root_hyper_sets[t];
    //         traverse_args[t]->leaf = NULL;
    //         traverse_args[t]->q_point = q[0];
    //         pthread_create(&thread_ids[t], NULL, hyper_binary_traverse, traverse_args[t]);
    //     }

    //     for (t = 0; t < MAX_THREADS; t++) {
    //         pthread_join(thread_ids[t], NULL);
    //     }

    //     leafs_size = 0;
    //     for (t = 0; t < MAX_THREADS; t++) {
    //         leafs_size += traverse_args[t]->leaf->set_size;
    //     }

    //     leafs = (struct hyper_set **)malloc(leafs_size * sizeof(struct hyper_set *));
    //     for (t = 0; t < MAX_THREADS; t++) {
    //         leafs[t] = traverse_args[t]->leaf;
    //     }

    // }
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

    leafs_size = 0;
    for (t = 0; t < MAX_THREADS; t++) {
        leafs_size += traverse_args[t]->leaf->set_size;
    }

    leafs = (int **)malloc(leafs_size * sizeof(int *));
    for (t = 0; t < MAX_THREADS; t++) {
        leafs[t] = traverse_args[t]->leaf->indexes;
    }
    // qsort(leafs, (size_t) leafs_size, d * sizeof(int), compare_rows);

    clock_gettime(CLOCK_MONOTONIC, &end);

    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Tree scanning finished in: %lf seconds!\n\n", elapsed);
    // printf("Tree scanning finished in: %lf seconds!\n\n", (double)(clock() - start_t) / (MAX_THREADS * CLOCKS_PER_SEC));

    // Print all leafs!
    // for (t = 0; t < MAX_THREADS; t++) {
    //     for (i = 0; i < traverse_args[t]->leaf->set_size; i++) {
    //         for (j = 0; j < traverse_args[t]->leaf->d; j++) {
    //             printf("%lf ", traverse_args[t]->leaf->neighbors[i][j]);
    //         }
    //         printf("\b\n");
    //     }
    //     printf("\n\n");
    // }

    free(traverse_args);
    free(c);
    free(q);
    free(root_hyper_sets);

    return 0;
}

int compare_rows(const void *a, const void *b) {
    return memcmp(a, b, 700 * sizeof(int));
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

    hyper_position = 0.0;
    for (j = 0; j < hyper_subset->d; j++) {
        hyper_position += hyper_subset->coeffs[j] * q_point[j];
    }

    if (hyper_position > hyper_subset->b_value) {
        return hyper_subset->f_ptr[0];
    } else {
        return hyper_subset->f_ptr[1];
    }
}

void *hyper_binary_split(void *hyper_subset_void) {

    int i, j, t;
    int init1, init2;
    double *random_point_1, *random_point_2;
    double *midpoint, *normal_vector;
    double beta, hyper_position;
    struct hyper_set *hyper_subset, *new_hyper_subset_1, *new_hyper_subset_2;
    pthread_t thread_ids[2];

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

    new_hyper_subset_1 = malloc(sizeof(struct hyper_set));
    new_hyper_subset_1->indexes = (int *)malloc(hyper_subset->set_size * sizeof(int));
    new_hyper_subset_1->neighbors = (double **)malloc(hyper_subset->set_size * sizeof(double *));

    new_hyper_subset_2 = malloc(sizeof(struct hyper_set));
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
        if (hyper_subset->set_size == 20000) {
            pthread_create(&thread_ids[0], NULL, hyper_binary_split, new_hyper_subset_1);
        } else {
            hyper_binary_split(new_hyper_subset_1);
        }
    } else {
        init1 = 0;
        new_hyper_subset_1->f_ptr = NULL;
        // Strict kNN search!
    }
    
    if (new_hyper_subset_2->set_size > hyper_subset->depth) {
        if (hyper_subset->set_size == 20000) {
            pthread_create(&thread_ids[1], NULL, hyper_binary_split, new_hyper_subset_2);
        } else {
            hyper_binary_split(new_hyper_subset_2);
        }
        init2 = 1;
    } else {
        init2 = 0;
        new_hyper_subset_2->f_ptr = NULL;
        // Strict kNN search!
    }

    if (init1 && hyper_subset->set_size == 20000) {
        pthread_join(thread_ids[0], NULL);
    }
    
    if (init2 && hyper_subset->set_size == 20000) {
        pthread_join(thread_ids[1], NULL);
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
