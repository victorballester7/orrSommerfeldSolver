#ifndef MISC_H
#define LINALG_H

#include <complex.h>
#include <openblas/lapacke.h>

void diagMatrixTimesMatrix(int d, double *diag, double *matrix, double *result);
void matrixmult(int d, double *A, double *B, double *C);
void matrixsum(int d, double *A, double *B, double *C);
void diff_4thorder(int d, double *f, double *df, double h);
void readUFromFile(double *y, double *U, char *filename);
void interpolate(int dimU, double *y, double *U, int d, double *y_interp,
                 double *U_interp);
int countLinesInFile(char *filename);
void printEVtoFile(int d, double complex *lambda, char *filenameEValues);
#endif // MISC_H
