#include <complex.h>
#include <openblas/lapack.h>
#include <openblas/lapacke.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../include/EVproblem.h"
#include "../include/misc.h"

void printVector(int d, double *v) {
  for (int i = 0; i < d; i++) {
    printf("%f ", v[i]);
  }
  printf("\n");
}

void printMatrix(int d, double *A) {
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      printf("%f ", A[i * d + j]);
    }
    printf("\n");
  }
}

void printcMatrix(int d, double complex *A) {
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      printf("%f + %fi ", creal(A[i * d + j]), cimag(A[i * d + j]));
    }
    printf("\n");
  }
}

double complex runUFromFile(int N, double Re, double complex var,
                            double complex beta, char *filenameU,
                            char *filenameEValues, int isTemporal) {
  double *y, *U, *Uy, *Uyy, *U_interp, *Uy_interp, *Uyy_interp;
  lapack_complex_double *L, *M;
  double *aux;

  // data for change of variables
  double y_max = 25.;
  double yl =
      5.; // position where we want half of the points to be distributed to
  double a = yl * y_max / (y_max - 2 * yl);
  double b = 1 + 2 * a / y_max;

  clock_t start, end;

  start = clock();
  double *y_cheb = (double *)malloc((N - 1) * sizeof(double));
  double *y_chebextended = (double *)malloc((N + 1) * sizeof(double));
  double *y_cheb_new = (double *)malloc((N - 1) * sizeof(double));
  double *dp = (double *)malloc((N + 1) * (N + 1) * 4 * sizeof(double));
  double *dp_clamped = (double *)malloc((N - 1) * (N - 1) * 4 * sizeof(double));
  double *dp_new = (double *)malloc((N - 1) * (N - 1) * 4 * sizeof(double));

  aux = (double *)malloc((N - 1) * sizeof(double));
  U_interp = (double *)malloc((N - 1) * sizeof(double));
  Uy_interp = (double *)malloc((N - 1) * sizeof(double));
  Uyy_interp = (double *)malloc((N - 1) * sizeof(double));

  L = (lapack_complex_double *)malloc((N - 1) * (N - 1) *
                                      sizeof(lapack_complex_double));
  M = (lapack_complex_double *)malloc((N - 1) * (N - 1) *
                                      sizeof(lapack_complex_double));

  // getChebMatrices(N, y_cheb, D1, D2, D3, D4);

  // applyClampedBC(N, y_cheb, D1, D2, D3, D4, y_cheb_clamped, D1_clamped,
  //                D2_clamped, D3_clamped, D4_clamped);

  getChebClamped(N, y_cheb, dp, dp_clamped);

  getChebMatrices(N, y_chebextended, dp, dp + (N + 1) * (N + 1),
                  dp + 2 * (N + 1) * (N + 1), dp + 3 * (N + 1) * (N + 1));

  // printf("D4\n");
  // printf("D4\n");
  // for (int k = 0; k< 4; k++){
  //   for (int i = 0; i < N + 1; i++) {
  //     for (int j = 0; j < N + 1; j++) {
  //       printf("%f ", dp[k * (N + 1) * (N + 1) + i * (N + 1) + j]);
  //     }
  //     printf("\n");
  //   }
  //   printf("\n");
  // }

  // for (int i = 0; i < (N - 1); i++) {
  //   for (int j = 0; j < N - 1; j++) {
  //     dp_clamped[i * (N - 1) + j] = dp[(i + 1) * (N + 1) + (j + 1)];
  //     dp_clamped[i * (N - 1) + j + (N - 1) * (N - 1)] =
  //         dp[(i + 1) * (N + 1) + (j + 1) + (N + 1) * (N + 1)];
  //   }
  // }

  // printf("D2_clamped\n");

  // printf("D4_clamped\n");
  // printMatrix(N - 1, dp_clamped + 3 * (N - 1) * (N - 1));

  changeOfVariable(N - 1, y_cheb, dp_clamped, dp_clamped + (N - 1) * (N - 1),
                   dp_clamped + 2 * (N - 1) * (N - 1),
                   dp_clamped + 3 * (N - 1) * (N - 1), y_cheb_new, dp_new,
                   dp_new + (N - 1) * (N - 1), dp_new + 2 * (N - 1) * (N - 1),
                   dp_new + 3 * (N - 1) * (N - 1), a, b);

  // changeOfVariable(N - 1, y_cheb_clamped, D1_clamped, D2_clamped, D3_clamped,
  //                  D4_clamped, y_cheb_new, D1_new, D2_new, D3_new, D4_new, a,
  //                  b);

  // printMatrix(N - 1, dp_new + 1 * (N - 1) * (N - 1));
  // printVector(N - 1, y_cheb_new);

  // read file
  int dimU = countLinesInFile(filenameU);
  // remove header lines
  dimU -= 3;

  y = (double *)malloc(dimU * sizeof(double));
  U = (double *)malloc(dimU * sizeof(double));
  readUFromFile(y, U, filenameU);

  Uy = (double *)malloc(dimU * sizeof(double));
  Uyy = (double *)malloc(dimU * sizeof(double));
  diff_4thorder(dimU, U, Uy, y[1] - y[0]);
  diff_4thorder(dimU, Uy, Uyy, y[1] - y[0]);

  // shifted chebyshev nodes
  // for (int i = 0; i < N - 1; i++) {
  //   aux[i] = y_cheb_new[i] - y_cheb_new[N - 2];
  // }

  interpolate(dimU, y, U, N - 1, y_cheb_new, U_interp);
  interpolate(dimU, y, Uy, N - 1, y_cheb_new, Uy_interp);
  interpolate(dimU, y, Uyy, N - 1, y_cheb_new, Uyy_interp);

  // printVector(N - 1, U_interp);
  // printVector(N - 1, Uy_interp);
  // printVector(N - 1, Uyy_interp);

  // lapack setup for the system L v = lambda M v
  // lambda = ALPHA / BETA
  lapack_complex_double *ALPHA, *BETA;
  // lapack_complex_double *EIGENVECTORS;
  double complex *lambda;
  ALPHA =
      (lapack_complex_double *)malloc((N - 1) * sizeof(lapack_complex_double));
  BETA =
      (lapack_complex_double *)malloc((N - 1) * sizeof(lapack_complex_double));
  // EIGENVECTORS = (lapack_complex_double *)malloc((N - 1) * (N - 1) *
  //
  // sizeof(lapack_complex_double));
  lambda = (double complex *)malloc((N - 1) * sizeof(double complex));

  if (isTemporal) {
    // printf("alpha = %f\n", var);
    // printf("beta = %f\n", beta);
    getOSMatricesTemporal(N - 1, Re, var, beta, dp_new + 1 * (N - 1) * (N - 1),
                          dp_new + 3 * (N - 1) * (N - 1), U_interp, Uyy_interp,
                          L, M);
    end = clock();

    printf("Time to get matrices: %f\n",
           (double)(end - start) / CLOCKS_PER_SEC);

    // int info =
    //     LAPACKE_zggev(LAPACK_ROW_MAJOR, 'N', 'V', N - 1, L, N - 1, M, N - 1,
    //                   ALPHA, BETA, NULL, N - 1, EIGENVECTORS, N - 1);
    int info = LAPACKE_zggev(LAPACK_ROW_MAJOR, 'N', 'N', N - 1, L, N - 1, M,
                             N - 1, ALPHA, BETA, NULL, N - 1, NULL, N - 1);

    start = clock();

    printf("Time to get eigenvalues: %f\n",
           (double)(start - end) / CLOCKS_PER_SEC);

    if (info > 0) {
      printf("Failed to compute eigenvalues.\n");
      return -1;
    }
    for (int i = 0; i < N - 1; i++) {
      lambda[i] = ALPHA[i] / BETA[i] / var;
    }
    printEVtoFile(N - 1, lambda, filenameEValues);

  } else {
    // getOSMatrices(N - 1, Re, var, beta, D2_new, D4_new, U_interp,
    // Uy_interp,
    // Uyy_interp, L, M);
  }

  // free memory
  // free(y_cheb);
  // free(D1);
  // free(D2);
  // free(D3);
  // free(D4);
  // free(y_cheb_new);
  // free(D1_new);
  // free(D2_new);
  // free(D3_new);
  // free(D4_new);
  // free(y_cheb_clamped);
  // free(D1_clamped);
  // free(D2_clamped);
  // free(D3_clamped);
  // free(D4_clamped);
  // free(aux);
  // free(y);
  // free(U);
  // free(U_interp);
  // free(Uy_interp);
  // free(Uyy_interp);
  // free(L);
  // free(M);

  return 0;
}

int main() {
  const int N = 200;
  const double Re = 800;
  const double complex alpha = 1. + 0.*I;
  const double complex beta = 0.;
  const double complex omega = 0.1;
  const double deltaStar = 1.7207876573;

  int isTemporal = 1;
  char filenameU[50] = "data/blasius.dat";
  char filenameEValues[50] = "data/EValues.dat";

  runUFromFile(N, Re / deltaStar, alpha / deltaStar, beta / deltaStar,
               filenameU, filenameEValues, isTemporal);

  return 0;
}
