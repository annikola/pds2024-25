#ifndef MAT_READ_WRITE_H
#define MAT_READ_WRITE_H

#include <stdlib.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

double **read_2D_array_from_matfile(const char *filename, const char *varname, size_t *c_size, size_t *d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);

#endif // MAT_READ_WRITE_H
