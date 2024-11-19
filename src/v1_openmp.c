#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "../include/knn_search.h"
#include "../include/mat_read_write.h"

#define MIN_ARGS 5
#define MAX_SPLIT_DEPTH 2

typedef struct {
    double *coeffs;
    double b_value;
    Point **all_points;
    Point **points;
    int set_size;
    int d;
    int knns;
    int depth;
    int delta;
    int split_depth;
} hyper_set;

int is_duplicate(double *point1, double *point2, int d);
void *hyper_binary_split(void *hyper_subset);

int main(int argc, char *argv[]) {

    int i, j;
    int knns, depth, delta;
    size_t c_size, d;
    double elapsed;
    double **c, **my_idx, **my_dst;
    clock_t start_t;
    struct timespec start, end;
    hyper_set *root_hyper_set;
    const char *filename, *varname;
    Point **c_points;
    
    if (argc < MIN_ARGS + 1) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    filename = argv[1];
    varname = argv[2];
    knns = atoi(argv[3]);
    depth = atoi(argv[4]);
    delta = atoi(argv[5]);

    if (depth > 100000 || delta > 3000) {
        printf("This will ruin you DON'T DO IT (Change me inside the source file if you insist...)!\n");
        return 0;
    }

    printf("Reading the corpus...\n");
    c = read_2D_array_from_matfile(filename, varname, &c_size, &d);

    c_points = (Point **)malloc(c_size * sizeof(Point *));
    _2D_array_to_points(c_points, c, c_size, d);
    for (i = 0; i < c_size; i++) {
        free(c[i]);
    }
    free(c);

    printf("Initializing the splitting process...\n");
    clock_gettime(CLOCK_MONOTONIC, &start);

    srand(time(0));

    root_hyper_set = malloc(sizeof(hyper_set));
    root_hyper_set->all_points = c_points;
    root_hyper_set->points = c_points;
    root_hyper_set->set_size = (int)c_size;
    root_hyper_set->d = (int)d;
    root_hyper_set->coeffs = NULL; // Technically a zero-vector...
    root_hyper_set->b_value = 0.0;
    root_hyper_set->knns = knns;
    root_hyper_set->depth = depth;
    root_hyper_set->delta = delta;
    root_hyper_set->split_depth = 0;
    hyper_binary_split(root_hyper_set);

    clock_gettime(CLOCK_MONOTONIC, &end);

    // Calculate the elapsed time in seconds
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("Multithreaded application finished in: %lf seconds!\n\n", elapsed);

    // printf("Fixing files format...\n");
    my_idx = (double **)malloc(c_size * sizeof(double *));
    my_dst = (double **)malloc(c_size * sizeof(double *));
    for (i = 0; i < c_size; i++) {
        my_idx[i] = (double *)malloc(knns * sizeof(double));
        my_dst[i] = (double *)malloc(knns * sizeof(double));
        for (j = 0; j < knns; j++) {
            my_idx[i][j] = (double)c_points[i]->neighbors[j].index;
            my_dst[i][j] = c_points[i]->neighbors[j].distance;
        }
    }

    printf("Writing to mat files...\n");
    // write_2D_array_to_matfile("datasets/my_c.mat", "C", c, c_size, d);
    write_2D_array_to_matfile("datasets/my_idx.mat", "iii", my_idx, c_size, knns);
    write_2D_array_to_matfile("datasets/my_dst.mat", "ddd", my_dst, c_size, knns);
    
    // free(c);
    free(my_idx);
    free(my_dst);
    free(root_hyper_set);

    return 0;
}

int is_duplicate(double *point1, double *point2, int d) {

    int i;

    for (i = 0; i < d; i++) {
        if (point1[i] != point2[i]) {
            return 0;
        }
    }

    return 1;
}

void *hyper_binary_split(void *hyper_subset_void) {

    int i, j, edge_points_size;
    double *C1, *C2, *C3;
    double *random_point_1, *random_point_2;
    double *midpoint, *normal_vector;
    double beta, hyper_position;
    hyper_set *hyper_subset, *new_hyper_subset_1, *new_hyper_subset_2;
    Point **edge_points;

    hyper_subset = (hyper_set *)hyper_subset_void;

    random_point_1 = hyper_subset->points[rand() % hyper_subset->set_size]->coordinates;
    random_point_2 = hyper_subset->points[rand() % hyper_subset->set_size]->coordinates;

    // RARE CASE!!! CHECK IT THOUGH!!!
    while (is_duplicate(random_point_1, random_point_2, hyper_subset->d)) {
        random_point_2 = hyper_subset->points[rand() % hyper_subset->set_size]->coordinates;
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

    new_hyper_subset_1 = malloc(sizeof(hyper_set));
    new_hyper_subset_1->points = (Point **)malloc(hyper_subset->set_size * sizeof(Point *));

    new_hyper_subset_2 = malloc(sizeof(hyper_set));
    new_hyper_subset_2->points = (Point **)malloc(hyper_subset->set_size * sizeof(Point *));

    new_hyper_subset_1->all_points = new_hyper_subset_2->all_points = hyper_subset->all_points;
    new_hyper_subset_1->set_size = new_hyper_subset_2->set_size = 0;
    new_hyper_subset_1->d = new_hyper_subset_2->d = hyper_subset->d;
    new_hyper_subset_1->knns = new_hyper_subset_2->knns = hyper_subset->knns;
    new_hyper_subset_1->depth = new_hyper_subset_2->depth = hyper_subset->depth;
    new_hyper_subset_1->delta = new_hyper_subset_2->delta = hyper_subset->delta;
    new_hyper_subset_1->split_depth = new_hyper_subset_2->split_depth = hyper_subset->split_depth + 1;
    edge_points = (Point **)malloc(hyper_subset->set_size * sizeof(Point *));
    edge_points_size = 0;
    for (i = 0; i < hyper_subset->set_size; i++) {
        hyper_position = 0.0;
        for (j = 0; j < hyper_subset->d; j++) {
            hyper_position += normal_vector[j] * hyper_subset->points[i]->coordinates[j];
        }

        if (hyper_position > beta) {
            new_hyper_subset_1->points[new_hyper_subset_1->set_size] = hyper_subset->points[i];
            new_hyper_subset_1->set_size++;
        }
        else {
            new_hyper_subset_2->points[new_hyper_subset_2->set_size] = hyper_subset->points[i];
            new_hyper_subset_2->set_size++;
        }

        if (fabs(hyper_position - beta) < hyper_subset->delta) {
            edge_points[edge_points_size] = hyper_subset->points[i];
            edge_points_size++;
        }
    }

    free(normal_vector);
    free(midpoint);

    if (new_hyper_subset_1->set_size < hyper_subset->knns || new_hyper_subset_2->set_size < hyper_subset->knns) {
        if (hyper_subset->set_size > hyper_subset->depth) {
            hyper_binary_split(hyper_subset_void);
        } else {
            C1 = (double *)malloc(hyper_subset->set_size * hyper_subset->d * sizeof(double));
            _points_to_2D_mono_array(C1, hyper_subset->points, hyper_subset->set_size, hyper_subset->d);
            knn_search(C1, C1, hyper_subset->set_size, hyper_subset->set_size, hyper_subset->d, hyper_subset->knns, hyper_subset->points, hyper_subset->points, hyper_subset->all_points, 0);
            free(C1);
        }
    } else {
        new_hyper_subset_1->points = realloc(new_hyper_subset_1->points, new_hyper_subset_1->set_size * sizeof(Point *));

        new_hyper_subset_2->points = realloc(new_hyper_subset_2->points, new_hyper_subset_2->set_size * sizeof(Point *));

        edge_points = realloc(edge_points, edge_points_size * sizeof(Point *));

        #pragma omp parallel sections
        {
            #pragma omp section
            {
                // printf("NHS1: %d\n", new_hyper_subset_1->set_size);
                if (new_hyper_subset_1->set_size > hyper_subset->depth) {
                    hyper_binary_split(new_hyper_subset_1);
                } else {
                    C1 = (double *)malloc(new_hyper_subset_1->set_size * new_hyper_subset_1->d * sizeof(double));
                    _points_to_2D_mono_array(C1, new_hyper_subset_1->points, new_hyper_subset_1->set_size, new_hyper_subset_1->d);
                    knn_search(C1, C1, new_hyper_subset_1->set_size, new_hyper_subset_1->set_size, new_hyper_subset_1->d, new_hyper_subset_1->knns, new_hyper_subset_1->points, new_hyper_subset_1->points, hyper_subset->all_points, 0);
                    free(C1);
                }
            }

            #pragma omp section
            {
                // printf("NHS2: %d\n", new_hyper_subset_2->set_size);
                if (new_hyper_subset_2->set_size > hyper_subset->depth) {
                    hyper_binary_split(new_hyper_subset_2);
                } else {
                    C2 = (double *)malloc(new_hyper_subset_2->set_size * new_hyper_subset_2->d * sizeof(double));
                    _points_to_2D_mono_array(C2, new_hyper_subset_2->points, new_hyper_subset_2->set_size, new_hyper_subset_2->d);
                    knn_search(C2, C2, new_hyper_subset_2->set_size, new_hyper_subset_2->set_size, new_hyper_subset_2->d, new_hyper_subset_2->knns, new_hyper_subset_2->points, new_hyper_subset_2->points, hyper_subset->all_points, 0);
                    free(C2);
                }
            }
        }

        if (edge_points_size >= hyper_subset->knns) {
            // printf("Stitching NH1 and NH2...\n");
            // printf("Number of edge points: %d\n", edge_points_size);
            C3 = (double *)malloc(edge_points_size * hyper_subset->d * sizeof(double));
            _points_to_2D_mono_array(C3, edge_points, edge_points_size, hyper_subset->d);
            knn_search(C3, C3, edge_points_size, edge_points_size, hyper_subset->d, hyper_subset->knns, edge_points, edge_points, hyper_subset->all_points, 1);
            free(C3);
        }
    }

    free(new_hyper_subset_1);
    free(new_hyper_subset_2);
    free(edge_points);
}
