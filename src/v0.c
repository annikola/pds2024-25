#include "../include/knn_search.h"
#include "../include/mat_read_write.h"

#define MIN_ARGS 5
#define MIN_FEATURE 1.0
#define MAX_FEATURE 100.0
#define MAX_GB 1

double random_double(double min, double max);

int main(int argc, char *argv[]) {

    int i, j;
    int c_size, q_size, d, knns;
    double elapsed;
    double *C, *Q;
    double **corpus, **query, **my_idx, **my_dst;
    struct timespec start, end;
    // const char *filename;
    Point **c_points, **q_points;
    
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

    // This is the format that the matlab_write accepts...
    corpus = (double **)malloc(c_size * sizeof(double *));
    for (i = 0; i < c_size; i++) {
        corpus[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            corpus[i][j] = random_double(MIN_FEATURE, MAX_FEATURE);
        }
    }

    query = (double **)malloc(q_size * sizeof(double *));
    for (i = 0; i < q_size; i++) {
        query[i] = (double *)malloc(d * sizeof(double));
        for (j = 0; j < d; j++) {
            query[i][j] = random_double(MIN_FEATURE, MAX_FEATURE);
        }
    }

    /* <-- THA FYGEI STO TELOS!!! */

    // printf("Reading the corpus...\n");
    // c = read_2D_array_from_matfile(filename, varname, &c_size, &d);

    c_points = (Point **)malloc(c_size * sizeof(Point *));
    _2D_array_to_points(c_points, corpus, c_size, d);
    for (i = 0; i < c_size; i++) {
        free(corpus[i]);
    }
    free(corpus);

    q_points = (Point **)malloc(q_size * sizeof(Point *));
    _2D_array_to_points(q_points, query, q_size, d);
    for (i = 0; i < q_size; i++) {
        free(query[i]);
    }
    free(query);

    C = (double *)malloc(c_size * d * sizeof(double));
    _points_to_2D_mono_array(C, c_points, c_size, d);
    Q = (double *)malloc(q_size * d * sizeof(double));
    _points_to_2D_mono_array(Q, q_points, q_size, d);

    printf("Initializing the k-NN search...\n");
    clock_gettime(CLOCK_MONOTONIC, &start);

    knn_search(C, Q, c_size, q_size, d, knns, c_points, q_points, NULL, 0);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("k-NN search finished in: %lf seconds!\n\n", elapsed);

    printf("Fixing files format...\n");
    my_idx = (double **)malloc(c_size * sizeof(double *));
    my_dst = (double **)malloc(c_size * sizeof(double *));
    for (i = 0; i < q_size; i++) {
        my_idx[i] = (double *)malloc(knns * sizeof(double));
        my_dst[i] = (double *)malloc(knns * sizeof(double));
        for (j = 0; j < knns; j++) {
            my_idx[i][j] = (double)q_points[i]->neighbors[j].index;
            my_dst[i][j] = q_points[i]->neighbors[j].distance;
        }
    }

    printf("Writing to mat files...\n");
    // write_2D_array_to_matfile("corpus.mat", "C", corpus, c_size, d);
    // write_2D_array_to_matfile("query.mat", "Q", query, q_size, d);
    write_2D_array_to_matfile("my_idx.mat", "iii", my_idx, q_size, knns);
    write_2D_array_to_matfile("my_dst.mat", "ddd", my_dst, q_size, knns);

    free(C);
    free(Q);
    free(my_idx);
    free(my_dst);

    return 0;
}

double random_double(double min, double max) {

    double scale;

    scale = rand() / (double) RAND_MAX;

    return min + scale * (max - min);
}
