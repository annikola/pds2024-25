#ifndef KNN_SEARCH_H
#define KNN_SEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <cblas.h>
#include <unistd.h>

#define Q_SPLIT 4
#define MAX_UPDATE_DEPTH 2
#define MAX_DEPTH 5000
#define MAX_PART 1000

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
    Point **all_points;
    Point **corpus_points;
    Point **query_part_points;
    int corpus_size;
    int query_part_size;
    int dimensions;
    int knns;
    int stitching;
} calculate_distances_args;

void _2D_array_to_points(Point **c_points, double **corpus, size_t corpus_size, size_t dimensions);
void _points_to_2D_mono_array(double *C, Point **c_points, size_t c_size, size_t d);
int qsort_compare(const void* a, const void* b);
int is_a_neighbor(Neighbor *neighbors, int knns, int current_index);
void update_neighbors(Point **all_points, Neighbor *neighbors, int dimensions, int knns, int index, int depth);
double *calculate_norms(double *matrix2D, int m_size, int d);
void *calculate_distances(void *args);
void knn_parallel_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **corpus_set_points, Point **query_set_points, Point **all_points, int stitching);
void knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **corpus_set_points, Point **query_set_points, Point **all_points, int stitching);

#endif // KNN_SEARCH_H
