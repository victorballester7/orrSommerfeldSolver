#include <complex.h>
#include <math.h>
#include <openblas/lapacke.h>
#include <stdio.h>
#include <stdlib.h>

#include "include/EVproblem.h"
#include "include/misc.h"

#define d1map2(i) (d1map[i] * d1map[i])
#define d1map3(i) (d1map[i] * d1map2(i))
#define d1map4(i) (d1map[i] * d1map3(i))

const double PI = 3.14159265358979323846;
void mapFromCheb(int d, double *y_cheb, double *y_real, double a, double b) {
  for (int i = 0; i < d; i++) {
    y_real[i] = a * (1. + y_cheb[i]) / (b - y_cheb[i]);
  }
}

void diffmapFromCheb(int d, double *y_cheb, double *y_real, double a, double b,
                     int order) {
  double tmp;
  if (order == 1) {
    for (int i = 0; i < d; i++) {
      tmp = a + y_cheb[i];
      y_real[i] = a * (b + 1.) / (tmp * tmp);
    }
  } else if (order == 2) {
    for (int i = 0; i < d; i++) {
      tmp = (a + y_cheb[i]);
      y_real[i] = -2. * a * (b + 1.) / (tmp * tmp * tmp);
    }
  } else if (order == 3) {
    for (int i = 0; i < d; i++) {
      tmp = (a + y_cheb[i]);
      y_real[i] = 6. * a * (b + 1.) / (tmp * tmp * tmp * tmp);
    }
  } else if (order == 4) {
    for (int i = 0; i < d; i++) {
      tmp = (a + y_cheb[i]);
      y_real[i] = -24. * a * (b + 1.) / (tmp * tmp * tmp * tmp * tmp);
    }
  }
}

// get chebysheb nodes with interval [-1,1] as well as the derivative matrix
// Note: the nodes are ordered from 1 to -1
void getCheb(int N, double *y, double *D) {
  for (int i = 0; i < N + 1; i++) {
    // naive approach
    // y[N - i] = cos(PI * i / N );

    // to get completely symmetric nodes around the origin due to floating point
    // errors. We use the identity cos(x) = sin(PI/2 - x)
    y[i] = sin(PI * (N - 2 * i) / (2 * N));
  }

  double ci, cj;
  int minusOnePow;

  for (int i = 0; i < N + 1; i++) {
    for (int j = 0; j < N + 1; j++) {
      if (i < N / 2.) {
        ci = (i == 0) ? 2.0 : 1.0;
        if (i == j) {
          if (i == 0) {
            D[i * (N + 1) + j] = (2.0 * N * N + 1.0) / 6.0;
          } else {
            // we change the denominator to avoid numerical instability for high
            // N 1 - y[i] * y[i] = sin^2(PI * i / N)
            D[i * (N + 1) + j] =
                -0.5 * y[i] / (sin(PI * i / N) * sin(PI * i / N));
          }
        } else {
          cj = (j == 0 || j == N) ? 2.0 : 1.0;
          minusOnePow = ((i + j) % 2 == 0) ? 1 : -1;

          // y[i] - y[j] = 2 * sin(PI * (-i + j) / (2 * N)) * sin(PI * (i + j) /
          // (2 * N))
          D[i * (N + 1) + j] =
              -0.5 * (ci / cj) * minusOnePow /
              (sin(PI * (i + j) / (2 * N)) * sin(PI * (i - j) / (2 * N)));
        }
      } else {
        // flipping technique (Donn & Solomonoff, 1995)
        D[i * (N + 1) + j] = -D[(N - i) * (N + 1) + (N - j)];
      }
    }
  }
}

void getChebMatrices(int N, double *y, double *D1, double *D2, double *D3,
                     double *D4) {
  getCheb(N, y, D1);
  matrixmult(N + 1, D1, D1, D2);
  matrixmult(N + 1, D2, D1, D3);
  matrixmult(N + 1, D3, D1, D4);
}

void printMatrix2(int d, double *A) {
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      printf("%f ", A[i * d + j]);
    }
    printf("\n");
  }
}

void printVector2(int d, double *v) {
  for (int i = 0; i < d; i++) {
    printf("%f ", v[i]);
  }
  printf("\n");
}

void applyClampedBC(int N, double *y_cheb, double *D1, double *D2, double *D3,
                    double *D4, double *y_cheb_new, double *D1_new,
                    double *D2_new, double *D3_new, double *D4_new) {
  // relation between the derivatives of p and h in the following formula:
  // p_i(x) = (1 - x^2)^2 / (1 - x_i^2)^2 h_i(x)
  // So we get:
  // p_i'(x)    = 1 / (1 - x_i^2)^2 [ (1 - x^2)^2 h_i'(x) - 4 x (1 - x^2) h_i(x)
  // ]

  // p_i''(x)   = 1 / (1 - x_i^2)^2 [ (1 - x^2)^2 h_i''(x) - 8 x (1 - x^2)
  // h_i'(x) + (12 x^2 - 4) h_i(x) ]

  // p_i'''(x)  = 1 / (1 - x_i^2)^2 [ (1 - x^2)^2 h_i'''(x) - 12 x (1 - x^2)
  // h_i''(x) + (36 x^2 - 12) h_i'(x) + 24 x h_i(x) ]

  // p_i''''(x) = 1 / (1 - x_i^2)^2 [ (1 - x^2)^2 h_i''''(x) - 16 x (1 - x^2)
  // h_i'''(x) + (72 x^2 - 24) h_i''(x) + 96 x h_i'(x) + 24 h_i(x) ]

  double tmp_i, tmp_j;
  // new matrices are of size (N + 1) - 2 = N - 1
  //   for (int i = 0; i < N - 1; i++) {
  //     y_cheb_new[i] = y_cheb[i + 1]; // interior chebyshev nodes

  //     // tmp1 = 1 / (1 - x_i^2)^2 = 1/ (sin(PI * i / N))^4
  //     tmp_i = sin(PI * (i + 1) / N);
  //     tmp_i *= tmp_i;
  //     tmp_i = 1. / (tmp_i * tmp_i);
  //     for (int j = 0; j < N - 1; j++) {
  //       // tmp_j = (1 - x^2) = (sin(PI * j / N))^2
  //       tmp_j = sin(PI * (j + 1) / N);
  //       tmp_j *= tmp_j;
  // #define D0ij (i == j ? 1. : 0.)
  // #define D1ij D1[(i + 1) * (N + 1) + j + 1]
  // #define D2ij D2[(i + 1) * (N + 1) + j + 1]
  // #define D3ij D3[(i + 1) * (N + 1) + j + 1]
  // #define D4ij D4[(i + 1) * (N + 1) + j + 1]
  // #define xj y_cheb[j + 1]
  //       D1_new[i * (N - 1) + j] = tmp_i * tmp_j * (tmp_j * D1ij - 4. * xj *
  //       D0ij); D2_new[i * (N - 1) + j] =
  //           tmp_i * (tmp_j * tmp_j * D2ij - 8. * xj * tmp_j * D1ij +
  //                    (12. * xj * xj - 4.) * D0ij);
  //       D3_new[i * (N - 1) + j] =
  //           tmp_i * (tmp_j * tmp_j * D3ij - 12. * xj * tmp_j * D2ij +
  //                    (36. * xj * xj - 12.) * D1ij + 24. * xj * D0ij);
  //       D4_new[i * (N - 1) + j] =
  //           tmp_i * (tmp_j * tmp_j * D4ij - 16. * xj * tmp_j * D3ij +
  //                    (72. * xj * xj - 24.) * D2ij + 96. * xj * D1ij + 24. *
  //                    D0ij);
  // #undef D0ij
  // #undef D1ij
  // #undef D2ij
  // #undef D3ij
  // #undef D4ij
  // #undef xj
  //     }
  //   }

  // using recursive relation from Weelfert (1997)
  // we want to find the derivatives of

  // D_ij^(k) = d^k/dx^k (alpha(x)/ alpha(x_j) * phi_j(x))|_{x = x_i}

  // where phi_j are the lagrange polynomials with chebysheb nodes and in our
  // case alpha(x) = 1 - x^2

  // off-diagonal elements satisfay (i != j)

  // D_ij^(k) = k / (x_i - x_j) * (c_i/c_j * D_ii^(k-1) - D_ij^(k-1))

  // where D_ij^(0) = delta_ij is the identity matrix
}

// Function to approximate alpha(m, x) - modify this as needed
double alpha(int m, double x) {
  if (m == 0) {
    return (1.0 - x * x) * (1.0 - x * x);
  } else if (m == 1) {
    return -4.0 * x * (1.0 - x * x);
  } else if (m == 2) {
    return -4. + 12. * x * x;
  } else if (m == 3) {
    return -24. * x;
  } else if (m == 4) {
    return -24.;
  } else {
    return 0.0;
  }
}

void getChebClamped(int N, double *y_cheb, double *dp, double *dp_clamp) {
  double *d = (double *)malloc((N - 1) * 4 * sizeof(double));
  double *d_clamp = (double *)malloc((N - 1) * 4 * sizeof(double));

  init(N, y_cheb, d_clamp);
  diag(N, d, d_clamp);
  offdiag(N, d, d_clamp, dp, dp_clamp);

  free(d);
}

void init(int N, double *x, double *d) {
  // Initialize d array
  // on the output d will contain the values of alpha^(m)(x_j)/alpha(x_j) for m
  // = 0, 1, 2, 3, 4 alpha(x) = (1 - x^2)^2 x contains the Chebyshev nodes xx
  // contains the evalutation 1 - x^2
  int dim = N - 1;
  double tmp;
  for (int i = 0; i < dim; i++) {
    x[i] = cos(PI * (i + 1) / N);
    tmp = sin(PI * (i + 1) / N) * sin(PI * (i + 1) / N);
    d[i] = -4.0 * x[i] / tmp;
    d[1 * dim + i] = (-4.0 + 12.0 * x[i] * x[i]) / (tmp * tmp);
    d[2 * dim + i] = 24.0 * x[i] / (tmp * tmp);
    d[3 * dim + i] = 24.0 / (tmp * tmp);
  }
}

// Computes diagonal elements of Dm
void diag(int N, double *d, double *d_clamp) {
  double tmp;
  int j;
  int dim = N - 1;

  for (int i = 0; i < dim; i++) {
    // for (int j = 0; j < N; j++) {
    j = (i + 1) % dim;
    // tmp = 1.0 / (x[i] - x[j]);
    d[i] = 0.0;
    d[dim + i] = 0.0;
    d[2 * dim + i] = 0.0;
    d[3 * dim + i] = 0.0;
    while (j != i) {
      tmp =
          2. * sin(PI * (j - i) / (2. * N)) * sin(PI * (i + j + 2) / (2. * N));
      tmp = 1.0 / tmp;

      // normal boundary conditions (with alpha(x) = 1)
      d[3 * dim + i] += 4.0 * tmp * d[2 * dim + i];
      d[2 * dim + i] += 3.0 * tmp * d[dim + i];
      d[dim + i] += 2.0 * tmp * d[i];
      d[i] += tmp;

      // clamped boundary conditions (with alpha(x) = (1 - x^2)^2)
      d_clamp[3 * dim + i] += 4.0 * tmp * d_clamp[2 * dim + i];
      d_clamp[2 * dim + i] += 3.0 * tmp * d_clamp[dim + i];
      d_clamp[dim + i] += 2.0 * tmp * d_clamp[i];
      d_clamp[i] += tmp;

      j = (j + 1) % dim;
    }
  }
}

// Computes off-diagonal elements of Dm
void offdiag(int N, double *d, double *d_clamp, double *dp, double *dp_clamp) {
  int dim = N - 1;
  int dimMat = dim * dim;
  double tmp;
  double ci, cj, ci_clamp, cj_clamp;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      if (i != j) {
        // it can be shown that the product c_i/c_j reduces (after a lot of
        // cancellations) to sin^2(PI * i / N) / sin^2(PI * j / N) * (-1)^(i -
        // j)
        ci_clamp = sin(PI * (i + 1) / N) * sin(PI * (i + 1) / N) *
                   ((i + 1) % 2 == 0 ? 1.0 : -1.0);
        cj_clamp = sin(PI * (j + 1) / N) * sin(PI * (j + 1) / N) *
                   ((j + 1) % 2 == 0 ? 1.0 : -1.0);
        ci_clamp = ci_clamp / cj_clamp;

        ci = ((i-j) % 2 == 0) ? 1.0 : -1.0;

        tmp = 2. * sin(PI * (j - i) / (2. * N)) *
              sin(PI * (i + j + 2) / (2. * N));
        tmp = 1.0 / tmp;

        // we do not store the identity matrix for D0

        // normal boundary conditions (with alpha(x) = 1)
        dp[i * dim + j] = tmp * ci;
        dp[i * dim + j + dimMat] = 2 * tmp * (ci * d[i] - dp[i * dim + j]);
        dp[i * dim + j + 2 * dimMat] =
            3 * tmp * (ci * d[i + dim] - dp[i * dim + j + 1 * dimMat]);
        dp[i * dim + j + 3 * dimMat] =
            4 * tmp * (ci * d[i + 2 * dim] - dp[i * dim + j + 2 * dimMat]);

        // clamped boundary conditions (with alpha(x) = (1 - x^2)^2)
        dp_clamp[i * dim + j] = tmp * ci_clamp;
        dp_clamp[i * dim + j + dimMat] =
            2 * tmp * (ci_clamp * d_clamp[i] - dp_clamp[i * dim + j]);
        dp_clamp[i * dim + j + 2 * dimMat] =
            3 * tmp *
            (ci_clamp * d_clamp[i + dim] - dp_clamp[i * dim + j + 1 * dimMat]);
        dp_clamp[i * dim + j + 3 * dimMat] =
            4 * tmp *
            (ci_clamp * d_clamp[i + 2 * dim] -
             dp_clamp[i * dim + j + 2 * dimMat]);
      } else {
        // normal boundary conditions (with alpha(x) = 1)
        dp[i * dim + j] = d[i];
        dp[i * dim + j + dimMat] = d[i + dim];
        dp[i * dim + j + 2 * dimMat] = d[i + 2 * dim];
        dp[i * dim + j + 3 * dimMat] = d[i + 3 * dim];

        // clamped boundary conditions (with alpha(x) = (1 - x^2)^2)
        dp_clamp[i * dim + j] = d_clamp[i];
        dp_clamp[i * dim + j + dimMat] = d_clamp[i + dim];
        dp_clamp[i * dim + j + 2 * dimMat] = d_clamp[i + 2 * dim];
        dp_clamp[i * dim + j + 3 * dimMat] = d_clamp[i + 3 * dim];
      }
    }
  }
}

void changeOfVariable(int d, double *y_cheb, double *D1, double *D2, double *D3,
                      double *D4, double *y_cheb_new, double *D1_new,
                      double *D2_new, double *D3_new, double *D4_new, double a,
                      double b) {
  mapFromCheb(d, y_cheb, y_cheb_new, a, b);

  double *d1map, *d2map, *d3map, *d4map;

  d1map = (double *)malloc(d * sizeof(double));
  d2map = (double *)malloc(d * sizeof(double));
  d3map = (double *)malloc(d * sizeof(double));
  d4map = (double *)malloc(d * sizeof(double));

  diffmapFromCheb(d, y_cheb_new, d1map, a, b, 1);
  diffmapFromCheb(d, y_cheb_new, d2map, a, b, 2);
  diffmapFromCheb(d, y_cheb_new, d3map, a, b, 3);
  diffmapFromCheb(d, y_cheb_new, d4map, a, b, 4);

  // printVector2(d, d1map);

  // printf("\n");
  // for (int i = 0; i < d; i++) {
  //   for (int j = 0; j < d; j++) {
  //     printf("%f ", D1[i * d + j]);
  //   }
  //   printf("\n");
  // }

  // printf("\n");

  // build of D1_new
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      D1_new[i * d + j] = d1map[i] * D1[i * d + j];
      // print/* f( */"%f ", D1_new[i * d + j]);
    }
    // printf("\n");
  }

  
  // printf("\n");


  // build of D2_new
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      D2_new[i * d + j] = d1map2(i) * D2[i * d + j] + d2map[i] * D1[i * d + j];
      // printf/* (" */%f ", D2_new[i * d + j]);
    }
    // printf("\n");
  }

  // build of D3_new
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      D3_new[i * d + j] = d1map3(i) * D3[i * d + j] +
                          3. * d1map[i] * d2map[i] * D2[i * d + j] +
                          d3map[i] * D1[i * d + j];
    }
  }

  // build of D4_new
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      D4_new[i * d + j] =
          d1map4(i) * D4[i * d + j] +
          6. * d1map2(i) * d2map[i] * D3[i * d + j] +
          (3. * d2map[i] * d2map[i] + 4. * d1map[i] * d3map[i]) *
              D2[i * d + j] +
          d4map[i] * D1[i * d + j];
    }
  }

  // free memory
  free(d1map);
  free(d2map);
  free(d3map);
  free(d4map);
}

void getOSMatricesTemporal(int d, double Re, double complex alpha,
                           double complex beta, double *D2, double *D4,
                           double *U, double *Uyy, lapack_complex_double *L,
                           lapack_complex_double *M) {
  double complex k2 = alpha * alpha + beta * beta;

  // print D4
  // printf("D4\n");
  // for (int i = 0; i < d; i++) {
  //   for (int j = 0; j < d; j++) {
  //     printf("%f ", D4[i * d + j]);
  //   }
  //   printf("\n");
  // }

  // printf("\n");

  // build of M
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      if (i == j) {
        M[i * d + j] = I * (k2 - D2[i * d + j]);
      } else {
        M[i * d + j] = -I * D2[i * d + j];
      }
      // printf("%f ", cimag(M[i * d + j]));
    }
    // printf("\n");
  }

  // printf("\n");
  // build of L
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      if (i == j) {
        L[i * d + j] =
            alpha * U[i] * M[i * d + j] + I * alpha * Uyy[i] +
            1. / Re * (D4[i * d + j] - 2 * k2 * D2[i * d + j] + k2 * k2);
      } else {
        L[i * d + j] = alpha * U[i] * M[i * d + j] +
                       1. / Re * (D4[i * d + j] - 2 * k2 * D2[i * d + j]);
      }
      // printf("%f r + %f i  ", creal(L[i * d + j]), cimag(L[i * d + j]));
    }
    // printf("\n");
  }
}

void solveEVproblemTemporal(double complex *L, double complex *M) {}
