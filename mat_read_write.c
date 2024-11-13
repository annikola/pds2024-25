#include <stdlib.h>
#include "/usr/local/MATLAB/R2024b/extern/include/mat.h"

double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d);
void write_2D_array_to_matfile(const char *filename, const char *array_name, double **_2D_array, int c_size, int d);

double **read_2D_array_from_matfile(const char *filename, size_t *c_size, size_t *d) {

    MATFile *pmat;
    mxArray *array_ptr;
    size_t i, j;
    double *data;
    double **c;
    const char *varname = "C"; // Name of the variable to read

    pmat = matOpen(filename, "r");
    if (pmat == NULL) {
        fprintf(stderr, "Error opening file data.mat\n");
        return NULL;
    }

    // Read a variable from the .mat file
    array_ptr = matGetVariable(pmat, varname);
    if (array_ptr == NULL) {
        fprintf(stderr, "Error reading variable %s from file\n", varname);
        matClose(pmat);
        return NULL;
    }

    // Check if the variable is of the expected type (for example, double matrix)
    if (mxIsDouble(array_ptr) && !mxIsComplex(array_ptr)) {
        // Get the data pointer and array dimensions
        data = mxGetPr(array_ptr);
        size_t rows = mxGetM(array_ptr);
        size_t cols = mxGetN(array_ptr);

        *c_size = rows;
        *d = cols;

        c = (double **)malloc(rows * sizeof(double *));
        for (i = 0; i < rows; i++) {
            c[i] = (double *)malloc(cols * sizeof(double));
        }

        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                c[i][j] = data[i + j * rows];
            }
        }
    } else {
        fprintf(stderr, "Variable %s is not a double matrix\n", varname);
        return NULL;
    }

    // Clean up
    mxDestroyArray(array_ptr); // Free the mxArray
    matClose(pmat);            // Close the MAT-file

    return c;
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
