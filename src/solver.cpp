#include "../include/solver.hpp"
#include "../libs/fastgl.hpp"
#include <cmath>
#include <complex>
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>

// Define complex number type
using complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

#define T1ijk (gaussWeights[k] * d2N[i][k] * d2N[j][k])
#define T2ijk (gaussWeights[k] * d2N[i][k] * N[j][k])
#define T3ijk (gaussWeights[k] * N[i][k] * N[j][k])
#define T4ijk (gaussWeights[k] * N[i][k] * d2U[k] * N[j][k])
#define T5ijk (gaussWeights[k] * U[k] * d2N[i][k] * N[j][k])
#define T6ijk (gaussWeights[k] * U[k] * N[i][k] * N[j][k])
#define T7ijk (gaussWeights[k] * N[i][k] * dN[j][k])
#define T8ijk (gaussWeights[k] * U[k] * dN[i][k] * d2N[j][k])
#define T9ijk (gaussWeights[k] * U[k] * N[i][k] * dN[j][k])

void CustomU::readFromFile(const std::string &filename,
                           std::vector<double> &xdata, uint colX,
                           std::vector<double> &ydata, uint colY,
                           uint numSkipHeaderLines) const {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  // Skip header lines
  std::string line;
  for (uint i = 0; i < numSkipHeaderLines; i++) {
    std::getline(file, line);
  }

  // Read data from file
  // add the data at colX to xdata and colY to data
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;

    while (iss >> token) {
      tokens.push_back(token);
    }

    // Check if we have enough columns
    if (tokens.size() > std::max(colX, colY)) {
      try {
        xdata.push_back(std::stod(tokens[colX]));
        ydata.push_back(std::stod(tokens[colY]));
      } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid number in file: " << line << std::endl;
      }
    } else if (!tokens.empty()) {
      std::cerr << "Not enough columns in line: " << line << std::endl;
    }
  }
  return;
}

double CustomU::interpolate(double z, const std::vector<double> &xdata,
                            const std::vector<double> &ydata) const {
  // Implement interpolation logic here
  // For simplicity, using linear interpolation
  uint index = 0;
  while (index < len_data - 1 && xdata[index] < z) {
    index++;
  }

  double t = (xdata[index] - z) / (xdata[index] - xdata[index - 1]);
  return ydata[index] - t * (ydata[index] - ydata[index - 1]);
}

std::vector<double> CustomU::diffData(const std::vector<double> &xdata,
                                      const std::vector<double> &ydata) const {
  std::vector<double> diff_data(ydata.size());
  // 1st order forward difference for all the elements except the last one:
  // double h = xdata[1] - xdata[0];
  size_t n = ydata.size();

  // 2nd order derivatives with non-uniform spacing
  if (n < 3) {
    throw std::invalid_argument("Input data must have at least 3 points.");
  }

  // boundary derivative
  double h0 = xdata[1] - xdata[0];
  double h1 = xdata[2] - xdata[1];
  diff_data[0] = -(2 * h0 + h1) / (h0 * (h0 + h1)) * ydata[0] +
                 (h0 + h1) / (h0 * h1) * ydata[1] -
                 h0 / (h1 * (h0 + h1)) * ydata[2];

  h0 = xdata[n - 1] - xdata[n - 2];
  h1 = xdata[n - 2] - xdata[n - 3];
  diff_data[n - 1] = (2 * h0 + h1) / (h0 * (h0 + h1)) * ydata[n - 1] -
                     (h0 + h1) / (h0 * h1) * ydata[n - 2] +
                     h0 / (h1 * (h0 + h1)) * ydata[n - 3];

  // interior derivatives
  for (uint i = 1; i < n - 1; i++) {
    h0 = xdata[i] - xdata[i - 1];
    h1 = xdata[i + 1] - xdata[i];
    diff_data[i] = -h1 / (h0 * (h0 + h1)) * ydata[i - 1] +
                   (h1 - h0) / (h0 * h1) * ydata[i] +
                   h0 / (h1 * (h0 + h1)) * ydata[i + 1];
  }

  // // 4th order centered difference for interior points
  // for (uint i = 2; i < n - 2; i++) {
  //   diff_data[i] =
  //       (ydata[i - 2] - 8 * ydata[i - 1] + 8 * ydata[i + 1] - ydata[i + 2]) /
  //       (12.0 * h);
  // }

  // // 4th order forward difference for the first element
  // diff_data[0] = (-25 * ydata[0] + 48 * ydata[1] - 36 * ydata[2] +
  //                 16 * ydata[3] - 3 * ydata[4]) /
  //                (12.0 * h);

  // // 4th order mixed difference for the second element
  // diff_data[1] = (-3 * ydata[0] - 10 * ydata[1] + 18 * ydata[2] - 6 *
  // ydata[3] +
  //                 ydata[4]) /
  //                (12.0 * h);

  // // 4th order mixed difference for the penultimate element
  // diff_data[n - 2] = (3 * ydata[n - 1] + 10 * ydata[n - 2] - 18 * ydata[n -
  // 3] +
  //                     6 * ydata[n - 4] - ydata[n - 5]) /
  //                    (12.0 * h);

  // // 4th order backward difference for the last element
  // diff_data[n - 1] =
  //     (25 * ydata[n - 1] - 48 * ydata[n - 2] + 36 * ydata[n - 3] -
  //      16 * ydata[n - 4] + 3 * ydata[n - 5]) /
  //     (12.0 * h);

  // for (uint i = 0; i < ydata.size() - 1; i++) {
  //   diff_data[i] = (ydata[i + 1] - ydata[i]) / (xdata[i + 1] - xdata[i]);
  // }

  return diff_data;
}

double OSSolver::getYPhysicalRegion(double y_standard) const {
  // Map the flow profile from the standard region [-1, 1] to the physical
  // region [a, b]
  return y_standard * (b - a) / 2.0 + (a + b) / 2.0;
}

void OSSolver::mapToStandardRegion(CustomU &Uprofile) {
  // Map the flow profile to the standard region [-1, 1]
  a = Uprofile.x_data[0];
  b = Uprofile.x_data[Uprofile.x_data.size() - 1];

  for (uint i = 0; i < Uprofile.x_data.size(); i++) {
    Uprofile.x_data[i] = 2.0 * (Uprofile.x_data[i] - a) / (b - a) - 1.0;
    // Uprofile.u_data[i] = 2.0 * (Uprofile.u_data[i] - a) / (b - a) - 1.0;
    // Uprofile.du_data[i] = 2.0 * (Uprofile.du_data[i] - a) / (b - a) - 1.0;
    // Uprofile.d2u_data[i] = 2.0 * (Uprofile.d2u_data[i] - a) / (b - a) - 1.0;
  }
}

// Evaluate Legendre polynomial of degree n at x using recursion formula
double OSSolver::getL(uint n, double x) const {
  if (n == 0)
    return 1.0;
  if (n == 1)
    return x;

  double p0 = 1.0;
  double p1 = x;
  double p2;

  for (uint i = 2; i <= n; i++) {
    p2 = ((2 * i - 1) * x * p1 - (i - 1) * p0) / i;
    p0 = p1;
    p1 = p2;
  }

  return p1;
}

// Internal shape function N_i(z) as defined in equation (10)
double OSSolver::getN(uint i, double z) const {
  assert(i >= 1 && i <= dimVS); // Ensure i is within valid range
  // Implementation based on equation (10) in the paper
  double factor = sqrt((2 * i + 3) / 2.0);

  // Alternative implementation using Corollary 2.5(a) from the paper
  double term1 =
      (getL(i + 3, z) - getL(i + 1, z)) / ((2 * i + 3) * (2 * i + 5));
  double term2 =
      (getL(i + 1, z) - getL(i - 1, z)) / ((2 * i + 1) * (2 * i + 3));

  return factor * (term1 - term2);
}

double OSSolver::getdN(uint i, double z) const {
  assert(i >= 1 && i <= dimVS); // Ensure i is within valid range
  double factor = sqrt((2 * i + 3) * 2.0);

  return 1.0 / factor * (getL(i + 2, z) - getL(i, z));
}

// Second derivative of shape function
double OSSolver::getd2N(uint i, double z) const {
  assert(i >= 1 && i <= dimVS); // Ensure i is within valid range
  // Implementation based on equation (c) in the paper after Definition 2.2
  double factor = sqrt((2 * i + 3) / 2.0);
  return factor * getL(i + 1, z);
}

void OSSolver::setGaussPointsWeights() {
  // Use Gauss-Legendre quadrature to get points and weights
  gaussPoints.resize(numQuadPoints);
  gaussWeights.resize(numQuadPoints);
  for (uint i = 1; i <= numQuadPoints; i++) {
    fastgl::QuadPair qp = fastgl::GLPair(numQuadPoints, i);
    gaussPoints[i - 1] = qp.x();
    gaussWeights[i - 1] = qp.weight;
  }
}

void OSSolver::setFunctions(Uprofile &Uprofile) {
  N.resize(dimVS, std::vector<double>(numQuadPoints));
  dN.resize(dimVS, std::vector<double>(numQuadPoints));
  d2N.resize(dimVS, std::vector<double>(numQuadPoints));
  U.resize(numQuadPoints);
  d2U.resize(numQuadPoints);
  double z;

  for (uint i = 0; i < dimVS; i++) {
    for (uint j = 0; j < numQuadPoints; j++) {
      z = gaussPoints[j];
      N[i][j] = getN(i + 1, z);
      dN[i][j] = getdN(i + 1, z) * jacobian;
      d2N[i][j] = getd2N(i + 1, z) * jacobian * jacobian;
      if (i == 0) {
        U[j] = Uprofile.getU(z);
        d2U[j] = Uprofile.getd2U(z);
      }
    }
  }
}

// Compute matrix entries using Gauss quadrature
complex OSSolver::Los(uint i, uint j) const {
  // Implementation based on equation (12) in the paper
  complex T1ij = 0.0, T2ij = 0.0, T3ij = 0.0;
  complex T4ij = 0.0, T5ij = 0.0, T6ij = 0.0;
  complex I = complex(0.0, 1.0);
  complex alpha = var;

  // Use Gauss quadrature for numerical integration
  for (uint k = 0; k < numQuadPoints; k++) {
    // T1: (D2N, D2N)
    T1ij += T1ijk;

    // T2: (D2N, N)
    T2ij += T2ijk;

    // T3: (N, N)
    T3ij += T3ijk;

    // T4: (U''N, N)
    T4ij += T4ijk;

    // T5: (U D2N, N)
    T5ij += T5ijk;

    // T6: (U N, N)
    T6ij += T6ijk;
  }

  // Combine terms according to equation (14)
  return (T1ij - 2.0 * k2 * T2ij + k2 * k2 * T3ij + I * alpha * re * T4ij -
          I * alpha * re * T5ij + I * alpha * k2 * re * T6ij) /
         jacobian;
}

complex OSSolver::M(uint i, uint j) const {
  // Implementation based on equation (15) in the paper
  complex T2ij = 0.0, T3ij = 0.0;
  complex I = complex(0.0, 1.0);

  // Use Gauss quadrature for numerical integration
  for (uint k = 0; k < numQuadPoints; k++) {
    // T2: (D2N, N)
    T2ij += T2ijk;

    // T3: (N, N)
    T3ij += T3ijk;
  }

  // we solve for c = omega / alpha (in order to bettwer control the line
  // centered at one for blasius profile) except when alpha is 0, then we solve
  // for omega
  return -1.0 * I * re * (T2ij - k2 * T3ij) / jacobian;
}

complex OSSolver::R0(uint i, uint j) const {
  // Implementation based on equation (16) in the paper
  complex T1ij = 0.0, T2ij = 0.0, T3ij = 0.0;
  complex I = complex(0.0, 1.0);
  complex beta2 = beta * beta;
  complex omega = var;

  for (uint k = 0; k < numQuadPoints; k++) {
    // T1: (D2N, D2N)
    T1ij += T1ijk;

    // T2: (D2N, N)
    T2ij += T2ijk;

    // T3: (N, N)
    T3ij += T3ijk;
  }

  return (-1. / re * T1ij + (2. / re * beta2 - I * omega) * T2ij +
          (I * omega * beta2 - 1. / re * beta2 * beta2) * T3ij) /
         jacobian;
}

complex OSSolver::R1(uint i, uint j) const {
  // Implementation based on equation (16) in the paper
  complex T4ij = 0.0, T5ij = 0.0, T6ij = 0.0;
  complex T7ij = 0.0, T8ij = 0.0;
  complex I = complex(0.0, 1.0);
  complex beta2 = beta * beta;
  complex omega = var;

  for (uint k = 0; k < numQuadPoints; k++) {
    // T4: (U''N, N)
    T4ij += T4ijk;

    // T5: (U D2N, N)
    T5ij += T5ijk;

    // T6: (U N, N)
    T6ij += T6ijk;

    // T7: (N, DN)
    T7ij += T7ijk;

    // T8: (U D2N, DN)
    T8ij += T8ijk;
  }

  return (-1. * I * T4ij + I * T5ij - I * beta2 * T6ij +
          (2. * I * omega - 4. * beta2 / re) * T7ij - 4. / re * T8ij) /
         jacobian;
}

complex OSSolver::R2(uint i, uint j) const {
  // Implementation based on equation (16) in the paper
  complex T2ij = 0.0, T9ij = 0.0;
  complex I = complex(0.0, 1.0);

  for (uint k = 0; k < numQuadPoints; k++) {
    // T2: (D2N, N)
    T2ij += T2ijk;

    // T9: (U N, DN)
    T9ij += T9ijk;
  }

  return (4. / re * T2ij + 2. * I * T9ij) / jacobian;
}

// Build the A and B matrices for the generalized eigenvalue problem
void OSSolver::buildMatricesTemporal() {
  A = Matrix::Zero(dimVS, dimVS);
  B = Matrix::Zero(dimVS, dimVS);

  for (uint i = 0; i < dimVS; i++) {
    for (uint j = 0; j < dimVS; j++) {
      A(i, j) = Los(i, j);
      B(i, j) = M(i, j);
    }
  }
}
// Build the A and B matrices for the generalized eigenvalue problem
void OSSolver::buildMatricesSpatial() {
  uint dimMat = 2 * dimVS;

  A = Matrix::Zero(dimMat, dimMat);
  B = Matrix::Zero(dimMat, dimMat);

  for (uint i = 0; i < dimVS; i++) {
    for (uint j = 0; j < dimVS; j++) {
      A(i, j) = R1(i, j);
      A(i, j + dimVS) = R0(i, j);
      B(i, j) = R2(i, j);

      if (i == j) {
        A(i + dimVS, j) = 1.0;
        B(i + dimVS, j + dimVS) = 1.0;
      }
    }
  }
}

Eigen::VectorXcd
OSSolver::computeEigenvector(const Eigen::VectorXcd &eigenvector_coeffs) const {
  Eigen::VectorXcd eigenvector(numQuadPoints);
  for (uint j = 0; j < numQuadPoints; j++) {
    eigenvector[j] = 0.0;
    for (uint i = 0; i < dimVS; i++) {
      eigenvector[j] += eigenvector_coeffs[i] * N[i][j];
    }
  }

  return eigenvector;
}

// Solve the generalized eigenvalue problem
Eigen::ComplexEigenSolver<Matrix> OSSolver::solve() const {
  // transpose 
  Eigen::MatrixXcd A_transpose = A.transpose();
  Eigen::MatrixXcd B_transpose = B.transpose();
  

  // invert matrix B
  Eigen::FullPivLU<Matrix> lu(B_transpose);
  Matrix B_inv = lu.inverse();
  // B_inv * A
  Matrix M = B_inv * A_transpose;
  // Eigenvalue decomposition
  Eigen::ComplexEigenSolver<Matrix> eig(M);
  if (eig.info() != 0) {
    throw std::runtime_error("Eigenvalue decomposition failed");
  }

  return eig;
}
