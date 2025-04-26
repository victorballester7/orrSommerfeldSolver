#ifndef EVPROBLEM_H
#define EVPROBLEM_H

#include <complex.h>
#include <openblas/lapacke.h>

void mapFromCheb(int d, double *y_cheb, double *y_real, double a, double b);

void diffmapFromCheb(int d, double *y_cheb, double *y_real, double a, double b,
                     int order);

// get chebysheb nodes with interval [-1,1] as well as the derivative matrix
void getCheb(int N, double *y, double *D);

void getChebMatrices(int N, double *y, double *D1, double *D2, double *D3,
                     double *D4);

void changeOfVariable(int d, double *y_cheb, double *D1, double *D2, double *D3,
                      double *D4, double *y_cheb_new, double *D1_new,
                      double *D2_new, double *D3_new, double *D4_new, double a,
                      double b);

void applyClampedBC(int d, double *y_cheb, double *D1, double *D2, double *D3,
               double *D4, double *y_cheb_new, double *D1_new, double *D2_new,
               double *D3_new, double *D4_new);

void getOSMatricesTemporal(int d, double Re, double complex alpha, double complex beta, double *D2, double *D4, double *U, double *Uyy, lapack_complex_double *L, lapack_complex_double *M);

// void solveEVproblemTemporal(double complex *L, double complex *M);
void init(int N, double *x, double *d);
void diag(int N, double *d, double *d_clamp);
void offdiag(int N, double *d, double *d_clamp, double *dp, double *dp_clamp);
void getChebClamped(int N, double *y_cheb, double *dp, double *dp_clamp);
#endif // EVPROBLEM_H
