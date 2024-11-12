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
    Point *query_part_points;
    int corpus_size;
    int query_part_size;
    int dimensions;
    int knns;
} calculate_distances_args;

typedef struct {
    double *coeffs;
    double b_value;
    Point *points;
    int set_size;
    int d;
    int depth;
} hyper_set;

void _2D_array_to_points(Point *c_points, double **corpus, size_t corpus_size, size_t dimensions);
void *hyper_binary_split(void *hyper_subset);
int is_duplicate(double *point1, double *point2, int d);
void knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns, double **my_idx, double **my_dst);
double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);

int main(int argc, char *argv[]) {

    int depth;
    size_t c_size, d;
    double elapsed;
    double **c, **q;
    clock_t start_t;
    struct timespec start, end;
    hyper_set *root_hyper_set;
    const char *filename;
    Point *c_points;

    /* THA FYGEI STO TELOS!!! --> */
    
    if (argc < MIN_ARGS) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    filename = argv[1];
    depth = atoi(argv[2]);

    printf("Reading the corpus...\n");
    c = read_2D_array_from_matfile("big_set.mat", &c_size, &d);
    q = c;

    // for (i = 0; i < c_size; i++) {
    //     for (j = 0; j < d; j++) {
    //         printf("%lf ", c[i][j]);
    //     }
    //     printf("\b\n");
    // }
    // printf("\n");

    /* <-- THA FYGEI STO TELOS!!! */

    printf("Initializing the splitting process...\n");
    clock_gettime(CLOCK_MONOTONIC, &start);

    srand(time(0));

    c_points = (Point *)malloc(c_size * sizeof(Point));
    _2D_array_to_points(c_points, c, c_size, d);

    root_hyper_set = malloc(sizeof(hyper_set));
    root_hyper_set->points = c_points;
    root_hyper_set->set_size = (int)c_size;
    root_hyper_set->d = (int)d;
    root_hyper_set->coeffs = NULL; // Technically a zero-vector...
    root_hyper_set->b_value = 0.0;
    root_hyper_set->depth = depth;
    hyper_binary_split(root_hyper_set);

    clock_gettime(CLOCK_MONOTONIC, &end);

    // Calculate the elapsed time in seconds
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("Multithreaded application finished in: %lf seconds!\n\n", elapsed);
    
    free(c);
    free(root_hyper_set);

    return 0;
}

void _2D_array_to_points(Point *c_points, double **c, size_t c_size, size_t d) {
    
    int i, j;

    for (i = 0; i < c_size; i++) {
        c_points[i].index = i + 1;
        c_points[i].coordinates = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            c_points[i].coordinates[j] = c[i][j];
        }
    }
}

void *hyper_binary_split(void *hyper_subset_void) {

    int i, j, t;
    double *random_point_1, *random_point_2;
    double *midpoint, *normal_vector;
    double beta, hyper_position;
    hyper_set *hyper_subset, *new_hyper_subset_1, *new_hyper_subset_2;
    pthread_t thread_ids[2];
    calculate_distances_args **cd_args;

    hyper_subset = (hyper_set *)hyper_subset_void;

    random_point_1 = hyper_subset->points[rand() % hyper_subset->set_size].coordinates;
    random_point_2 = hyper_subset->points[rand() % hyper_subset->set_size].coordinates;

    // RARE CASE!!! CHECK IT THOUGH!!!
    // while (is_duplicate(random_point_1, random_point_2, hyper_subset->d)) {
    //     random_point_2 = hyper_subset->neighbors[rand() % hyper_subset->set_size];
    // }

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

    new_hyper_subset_1 = malloc(sizeof(hyper_set));
    new_hyper_subset_1->points = (Point *)malloc(hyper_subset->set_size * sizeof(Point));

    new_hyper_subset_2 = malloc(sizeof(hyper_set));
    new_hyper_subset_2->points = (Point *)malloc(hyper_subset->set_size * sizeof(Point));

    new_hyper_subset_1->depth = new_hyper_subset_2->depth = hyper_subset->depth;
    new_hyper_subset_1->d = new_hyper_subset_2->d = hyper_subset->d;
    new_hyper_subset_1->set_size = new_hyper_subset_2->set_size = 0;
    for (i = 0; i < hyper_subset->set_size; i++) {
        hyper_position = 0.0;
        for (j = 0; j < hyper_subset->d; j++) {
            hyper_position += normal_vector[j] * hyper_subset->points[i].coordinates[j];
        }

        if (hyper_position > beta) {
            new_hyper_subset_1->points[new_hyper_subset_1->set_size] = hyper_subset->points[i];
            new_hyper_subset_1->set_size++;
        }
        else {
            new_hyper_subset_2->points[new_hyper_subset_2->set_size] = hyper_subset->points[i];
            new_hyper_subset_2->set_size++;
        }
    }

    new_hyper_subset_1->points = realloc(new_hyper_subset_1->points, new_hyper_subset_1->set_size * sizeof(Point));

    new_hyper_subset_2->points = realloc(new_hyper_subset_2->points, new_hyper_subset_2->set_size * sizeof(Point));

    printf("NHS1: %d\n", new_hyper_subset_1->set_size);
    if (new_hyper_subset_1->set_size > hyper_subset->depth) {
        hyper_binary_split(new_hyper_subset_1);
    } else {
        ;
        // knn_search(C, Q, c_size, q_size, d, knns, my_idx, my_dst);
    }
    
    printf("NHS2: %d\n", new_hyper_subset_2->set_size);
    if (new_hyper_subset_2->set_size > hyper_subset->depth) {
        hyper_binary_split(new_hyper_subset_2);
    } else {
        ;
        // Strict kNN search!
    }

    // Edw exw ypologisei hdh dyo strict knns...
    // Apo edw kai katw kanw to stitch...

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
