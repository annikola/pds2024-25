#include <stdio.h>
#include <cblas.h>

int main() {
    int n = 3;
    double a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9}; // 3x3 matrix
    double x[3] = {1, 1, 1};                   // 3x1 vector
    double y[3] = {0, 0, 0};                   // Result vector

    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, a, n, x, 1, 0.0, y, 1);

    printf("Result: ");
    for (int i = 0; i < n; i++) {
        printf("%f ", y[i]);
    }
    printf("\n");

    return 0;
}
