#include <stdio.h>
#include <stdlib.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"
#include "../include/mat_read_write.h"

#define MIN_ARGS 4

int main(int argc, char *argv[]) {

    int i, j, neighbors_found;
    int *hash_table;
    size_t my_idx_size, my_idx_d, true_idx_size, true_idx_d;
    double **my_idx, **true_idx;
    const char *filename1, *filename2, *varname1, *varname2;
    
    if (argc < MIN_ARGS + 1) {
        printf("Not enough arguments provided!\n");
        return 0;
    }

    filename1 = argv[1];
    varname1 = argv[2];
    filename2 = argv[3];
    varname2 = argv[4];

    my_idx = read_2D_array_from_matfile(filename1, varname1, &my_idx_size, &my_idx_d);
    true_idx = read_2D_array_from_matfile(filename2, varname2, &true_idx_size, &true_idx_d);

    if (my_idx_size != true_idx_size || my_idx_d != true_idx_d) {
        printf("Wrong dimensions for your idx file or the reference file!\n");
        return 1;
    }

    neighbors_found = 0;
    for (i = 0; i < my_idx_size; i++) {

        hash_table = (int *)calloc(my_idx_size, sizeof(int));
        for (j = 0; j < my_idx_d; j++) {
            hash_table[(int)my_idx[i][j] % my_idx_size] = 1;
        }

        for (j = 0; j < my_idx_d; j++) {
            if (hash_table[(int)true_idx[i][j] % my_idx_size] == 1) {
                neighbors_found++;
            }
        }

        free(hash_table);
    }

    if (!neighbors_found) {
        printf("Total failure!\n");
    } else {
        printf("Recall percentage: %lf%%\n", (double)100 * neighbors_found / (my_idx_size * my_idx_d));
    }

    free(my_idx);
    free(true_idx);

    return 0;
}
