#include "include/misc.h"

#include <openblas/lapack.h>
#include <stdio.h>
#include <stdlib.h>

void diagMatrixTimesMatrix(int d, double *diag, double *matrix,
                           double *result) {
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      result[i * d + j] = diag[i] * matrix[i * d + j];
    }
  }
}

void matrixmult(int d, double *A, double *B, double *C) {
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      C[i * d + j] = 0;
      for (int k = 0; k < d; k++) {
        C[i * d + j] += A[i * d + k] * B[k * d + j];
      }
    }
  }
}

void matrixsum(int d, double *A, double *B, double *C) {
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      C[i * d + j] = A[i * d + j] + B[i * d + j];
    }
  }
}

void diff_4thorder(int d, double *f, double *df, double h) {
  for (int i = 0; i < d; i++) {
    if (i == 0) {
      df[i] = (-25. * f[i] + 48. * f[i + 1] - 36. * f[i + 2] + 16. * f[i + 3] -
               3. * f[i + 4]) /
              (12. * h);
    } else if (i == 1) {
      df[i] = (-3. * f[i - 1] - 10. * f[i] + 18. * f[i + 1] - 6. * f[i + 2] +
               f[i + 3]) /
              (12. * h);
    } else if (i == d - 2) {
      df[i] = (3. * f[i + 1] + 10. * f[i] - 18. * f[i - 1] + 6. * f[i - 2] -
               f[i - 3]) /
              (12. * h);
    } else if (i == d - 1) {
      df[i] = (25. * f[i] - 48. * f[i - 1] + 36. * f[i - 2] - 16. * f[i - 3] +
               3. * f[i - 4]) /
              (12. * h);
    } else {
      df[i] = (f[i - 2] - 8. * f[i - 1] + 8. * f[i + 1] - f[i + 2]) / (12. * h);
    }
  }
}

int countLinesInFile(char *filename) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error opening file");
    return -1;
  }

  int lines = 0;
  char ch;
  while ((ch = fgetc(file)) != EOF) {
    if (ch == '\n') {
      lines++;
    }
  }

  fclose(file);
  return lines;
}

void readUFromFile(double *y, double *U, char *filename) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error: file not found\n");
    exit(1);
  }

  // skip first 3 lines, then scan something of the form %lf %lf %lf **** \n
  // till the end of the file, which correspond to x, y, U. You can ignore the
  // other columns
  char buffer[1024];
  for (int i = 0; i < 3; i++) {
    if (!fgets(buffer, sizeof(buffer), file)) {
      fclose(file);
      return;
    }
  }

  int dimU = 0;
  double x, y_val, U_val;
  while (fgets(buffer, sizeof(buffer), file)) {
    if (sscanf(buffer, "%lf %lf %lf", &x, &y_val, &U_val) == 3) {
      y[dimU] = y_val;
      U[dimU] = U_val;
      dimU++;
    }
  }

  fclose(file);

  return;
}

void interpolate(int dimU, double *y, double *U, int d, double *y_interp,
                 double *U_interp) {
  for (int i = 0; i < d; i++) {
    double yi = y_interp[i];

    // If yi is out of bounds, extrapolate using the nearest point
    if (yi <= y[0]) {
      U_interp[i] = U[0];
    } else if (yi >= y[dimU - 1]) {
      U_interp[i] = U[dimU - 1];
    } else {
      // Find the interval [y[j], y[j+1]] that contains yi
      int j = 0;
      while (j < dimU - 1 && y[j + 1] < yi) {
        j++;
      }

      // Perform linear interpolation: U_interp = U1 + (U2 - U1) * (yi - y1) /
      // (y2 - y1)
      double y1 = y[j], y2 = y[j + 1];
      double U1 = U[j], U2 = U[j + 1];
      U_interp[i] = U1 + (U2 - U1) * (yi - y1) / (y2 - y1);
    }
  }
}

void printEVtoFile(int d, double complex *lambda, char *filenameEValues) {
  FILE *fileEValues = fopen(filenameEValues, "w");
  // FILE *fileEVectors = fopen(filenameEVectors, "w");

  if (fileEValues == NULL) {
    printf("Error: file not found\n");
    exit(1);
  }

  for (int i = 0; i < d; i++) {
    fprintf(fileEValues, "%f %f\n", creal(lambda[i]), cimag(lambda[i]));
  }

  // for (int j = 0; j < d; j++) {
  //   for (int i = 0; i < d; i++) {
  //     fprintf(fileEVectors, "%f %f ", creal(EVs[i * d + j]),
  //             cimag(EVs[i * d + j]));
  //   }
  //   fprintf(fileEVectors, "\n");
  // }

  fclose(fileEValues);
  // fclose(fileEVectors);
}
