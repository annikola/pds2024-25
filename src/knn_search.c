#include "../include/knn_search.h"

void _2D_array_to_points(Point **c_points, double **c, size_t c_size, size_t d) {
    
    int i, j;

    for (i = 0; i < c_size; i++) {
        c_points[i] = (Point *)malloc(sizeof(Point));
        c_points[i]->coordinates = (double *)malloc(d * sizeof(double));
        c_points[i]->index = i + 1; // Ti kai an mpei stin apo panw grammi?
        for (j = 0; j < d; j++) {
            c_points[i]->coordinates[j] = c[i][j];
        }
    }
}

void _points_to_2D_mono_array(double *C, Point **c_points, size_t c_size, size_t d) {
    
    int i, j;

    for (i = 0; i < c_size; i++) {
        for (j = 0; j < d; j++) {
            C[i * d + j] = c_points[i]->coordinates[j];
        }
    }
}

int qsort_compare(const void *a, const void *b) {
    const Neighbor *neighbor1 = (const Neighbor *)a;
    const Neighbor *neighbor2 = (const Neighbor *)b;

    if (neighbor1->distance < neighbor2->distance) {
        return -1;
    } else if (neighbor1->distance > neighbor2->distance) {
        return 1;
    } else {
        return 0;
    }
}

int is_a_neighbor(Neighbor *neighbors, int knns, int current_index) {
    
    int i;

    for (i = 0; i < knns; i++) {
        if (neighbors[i].index == current_index) {
            return 1;
        }
    }

    return 0;
}

void update_neighbors(Point **all_points, Neighbor *neighbors, int dimensions, int knns, int index, int depth) {

    int i, j;
    double distance;

    if (depth == MAX_UPDATE_DEPTH) {
        return;
    }

    for (i = 0; i < knns; i++) {
        distance = 0.0;
        for (j = 0; j < dimensions; j++) {
            distance += pow(all_points[neighbors[i].index - 1]->coordinates[j] - all_points[index - 1]->coordinates[j], 2);
        }
        distance = sqrt(distance);
        if (!is_a_neighbor(all_points[neighbors[i].index - 1]->neighbors, knns, index) && (distance < all_points[neighbors[i].index - 1]->neighbors[knns - 1].distance)) {
            all_points[neighbors[i].index - 1]->neighbors[knns - 1].index = index;
            all_points[neighbors[i].index - 1]->neighbors[knns - 1].distance = distance;
            qsort(all_points[neighbors[i].index - 1]->neighbors, (size_t)knns, sizeof(Neighbor), qsort_compare);
        }
    }

    for (i = 0; i < knns; i++) {
        update_neighbors(all_points, all_points[neighbors[i].index - 1]->neighbors, dimensions, knns, index, depth + 1);
    }
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

    int i, j;
    double curr_distance;
    double *c_norms, *q_norms, *matrixmul;
    calculate_distances_args *cd_args;
    Neighbor *curr_neighbors;

    cd_args = (calculate_distances_args *)args;

    matrixmul = (double *)malloc(cd_args->corpus_size * cd_args->query_part_size * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, cd_args->query_part_size, cd_args->corpus_size, cd_args->dimensions, -2.0, cd_args->query_part, cd_args->dimensions, cd_args->corpus, cd_args->dimensions, 0.0, matrixmul, cd_args->corpus_size);

    c_norms = calculate_norms(cd_args->corpus, cd_args->corpus_size, cd_args->dimensions);
    q_norms = calculate_norms(cd_args->query_part, cd_args->query_part_size, cd_args->dimensions);

    if (!cd_args->stitching) {
        for (i = 0; i < cd_args->query_part_size; i++) {
            cd_args->query_part_points[i]->neighbors = (Neighbor *)malloc(cd_args->corpus_size * sizeof(Neighbor));
            for (j = 0; j < cd_args->corpus_size; j++) {
                cd_args->query_part_points[i]->neighbors[j].index = cd_args->corpus_points[j]->index; // HMM...

                curr_distance = q_norms[i] + c_norms[j] + matrixmul[i * cd_args->corpus_size + j];
                if (curr_distance < 0) {
                    curr_distance = 0.0;
                } else {
                    curr_distance = sqrt(curr_distance);
                }
                cd_args->query_part_points[i]->neighbors[j].distance = curr_distance;
            }
        }

        for (i = 0; i < cd_args->query_part_size; i++) {
            qsort(cd_args->query_part_points[i]->neighbors, (size_t)cd_args->corpus_size, sizeof(Neighbor), qsort_compare);
            curr_neighbors = cd_args->query_part_points[i]->neighbors;
            cd_args->query_part_points[i]->neighbors = (Neighbor *)malloc(cd_args->knns * sizeof(Neighbor));
            memcpy(cd_args->query_part_points[i]->neighbors, curr_neighbors, cd_args->knns * sizeof(Neighbor));
            free(curr_neighbors);
        }
    } else {
        for (i = 0; i < cd_args->query_part_size; i++) {
            for (j = 0; j < cd_args->corpus_size; j++) {
                curr_distance = q_norms[i] + c_norms[j] + matrixmul[i * cd_args->corpus_size + j];
                if (curr_distance < 0) {
                    curr_distance = 0.0;
                } else {
                    curr_distance = sqrt(curr_distance);
                }

                if (!is_a_neighbor(cd_args->query_part_points[i]->neighbors, cd_args->knns, cd_args->corpus_points[j]->index) && (curr_distance < cd_args->query_part_points[i]->neighbors[cd_args->knns - 1].distance)) {
                    cd_args->query_part_points[i]->neighbors[cd_args->knns - 1].index = cd_args->corpus_points[j]->index;
                    cd_args->query_part_points[i]->neighbors[cd_args->knns - 1].distance = curr_distance;
                    qsort(cd_args->query_part_points[i]->neighbors, (size_t)cd_args->knns, sizeof(Neighbor), qsort_compare);
                    update_neighbors(cd_args->all_points, cd_args->query_part_points[i]->neighbors, cd_args->dimensions, cd_args->knns, cd_args->corpus_points[j]->index, 0);
                }
            }
        }
    }

    free(matrixmul);
    free(c_norms);
    free(q_norms);
}

void knn_parallel_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **corpus_set_points, Point **query_set_points, Point **all_points, int stitching) {
    
    int i, j, q_part_size, split_factor;
    double **q_parts;
    pthread_t q_thread_ids[Q_SPLIT];
    calculate_distances_args **cd_args;
    Point ***set_points_parts;

    // THEORITIKA DEN EINAI TELEIO GIA C != Q
    if (c_size <= 1 && q_size <= 1) {
        return;
    }
    
    q_part_size = (int)q_size / Q_SPLIT;
    if (q_part_size < 1) {
        q_part_size = 1;
        split_factor = q_size;
    } else {
        split_factor = Q_SPLIT;
    }

    set_points_parts = (Point ***)malloc(split_factor * sizeof(Point **));
    q_parts = (double **)malloc(split_factor * sizeof(double *));
    for (i = 0; i < split_factor; i++) {
        set_points_parts[i] = query_set_points + i * q_part_size;
        q_parts[i] = Q + i * d * q_part_size;
    }

    cd_args = (calculate_distances_args **)malloc(split_factor * sizeof(calculate_distances_args *));
    for (i = 0; i < split_factor; i++) {
        cd_args[i] = malloc(sizeof(calculate_distances_args));
        cd_args[i]->corpus = C;
        cd_args[i]->query_part = q_parts[i];
        cd_args[i]->corpus_size = c_size;
        if (i == split_factor - 1 && q_size > Q_SPLIT) {
            cd_args[i]->query_part_size = q_part_size + q_size % Q_SPLIT; // Sto teleutaio part prostithetai kai to ypoloipo...
        } else {
            cd_args[i]->query_part_size = q_part_size;
        }
        cd_args[i]->corpus_points = corpus_set_points;
        cd_args[i]->all_points = all_points;
        cd_args[i]->query_part_points = set_points_parts[i];
        cd_args[i]->dimensions = d;
        cd_args[i]->knns = knns;
        cd_args[i]->stitching = stitching;
        pthread_create(&q_thread_ids[i], NULL, calculate_distances, cd_args[i]);
    }

    for (i = 0; i < split_factor; i++) {
        pthread_join(q_thread_ids[i], NULL);
    }

    free(set_points_parts);
    free(q_parts);
    for (i = 0; i < split_factor; i++) {
        free(cd_args[i]);
    }
    free(cd_args);
}

void knn_search(double *C, double *Q, int c_size, int q_size, int d, int knns, Point **corpus_set_points, Point **query_set_points, Point **all_points, int stitching) {

    int i, corpus_size, query_size, parallel_parts;
    double **q_parts;
    Point ***set_points_parts;

    if (q_size <= MAX_DEPTH) {
        knn_parallel_search(C, Q, c_size, q_size, d, knns, corpus_set_points, query_set_points, all_points, stitching);
        return;
    }

    printf("Parting large set!\n");
    if (q_size % MAX_PART < knns) {
        parallel_parts = q_size / MAX_PART;
    }
    else {
        parallel_parts = q_size / MAX_PART + 1;
    }

    set_points_parts = (Point ***)malloc(parallel_parts * sizeof(Point **));
    q_parts = (double **)malloc(parallel_parts * sizeof(double *));
    for (i = 0; i < parallel_parts; i++) {
        set_points_parts[i] = query_set_points + i * MAX_PART;
        q_parts[i] = Q + i * d * MAX_PART;
    }

    printf("k-NN search started!\n");
    for (i = 0; i < parallel_parts; i++) {
        if (i == parallel_parts - 1) {
            if (q_size % MAX_PART < knns) {
                knn_parallel_search(C, q_parts[i], c_size, MAX_PART + q_size % MAX_PART, d, knns, corpus_set_points, set_points_parts[i], all_points, stitching);
            } else {
                knn_parallel_search(C, q_parts[i], c_size, q_size % MAX_PART, d, knns, corpus_set_points, set_points_parts[i], all_points, stitching);
            }
        } else {
            knn_parallel_search(C, q_parts[i], c_size, MAX_PART, d, knns, corpus_set_points, set_points_parts[i], all_points, stitching);
        }   
    }
    printf("k-NN search finished!\n");

    free(q_parts);
    free(set_points_parts);
}
