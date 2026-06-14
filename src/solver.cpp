#include "../include/solver.hpp"
#include "../libs/fastgl.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <limits>
#include <stdexcept>

// Define complex number type
using complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;
using RealMatrix = Eigen::MatrixXd;
using RealVector = Eigen::VectorXd;

extern "C" {
void zggev_(char *jobvl, char *jobvr, int *n, std::complex<double> *A, int *lda,
            std::complex<double> *B, int *ldb, std::complex<double> *alpha,
            std::complex<double> *beta, std::complex<double> *VL, int *ldvl,
            std::complex<double> *VR, int *ldvr, std::complex<double> *work,
            int *lwork, double *rwork, int *info);
}

#define T1ijk (factor * d2N[i][k] * d2N[j][k])
#define T2ijk (factor * d2N[i][k] * N[j][k])
#define T3ijk (factor * N[i][k] * N[j][k])
#define T4ijk (factor * N[i][k] * d2U[k] * N[j][k])
#define T5ijk (factor * U[k] * d2N[i][k] * N[j][k])
#define T6ijk (factor * U[k] * N[i][k] * N[j][k])
#define T7ijk (factor * dN[i][k] * N[j][k])
#define T8ijk (factor * U[k] * d2N[i][k] * dN[j][k])
#define T9ijk (factor * U[k] * dN[i][k] * N[j][k])

void CustomU::readFromFile(const std::string &filename,
                           std::vector<double> &xdata, uint colX,
                           std::vector<double> &ydata, uint colY,
                           uint numSkipHeaderLines) const {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Error opening file: " + filename);
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
  if (xdata.size() < 2 || ydata.size() < 2 || xdata.size() != ydata.size()) {
    throw std::invalid_argument(
        "Interpolation requires at least two aligned data points.");
  }
  if (z <= xdata.front()) {
    return ydata.front();
  }
  if (z >= xdata.back()) {
    return ydata.back();
  }

  auto upper = std::lower_bound(xdata.begin(), xdata.end(), z);
  size_t index = static_cast<size_t>(std::distance(xdata.begin(), upper));
  size_t index0 = index - 1;

  double x0 = xdata[index0];
  double x1 = xdata[index];
  double y0 = ydata[index0];
  double y1 = ydata[index];
  double t = (z - x0) / (x1 - x0);
  return y0 + t * (y1 - y0);
}

std::vector<double> CustomU::buildNaturalSplineSecondDerivatives(
    const std::vector<double> &xdata, const std::vector<double> &ydata) const {
  const size_t n = xdata.size();
  if (n < 3 || ydata.size() != n) {
    throw std::invalid_argument(
        "Spline construction requires at least 3 aligned data points.");
  }

  std::vector<double> y2(n, 0.0);
  std::vector<double> u(n, 0.0);

  for (size_t i = 1; i < n - 1; ++i) {
    const double hPrev = xdata[i] - xdata[i - 1];
    const double hNext = xdata[i + 1] - xdata[i];
    const double hTot = xdata[i + 1] - xdata[i - 1];
    if (hPrev <= 0.0 || hNext <= 0.0 || hTot <= 0.0) {
      throw std::invalid_argument(
          "Spline construction requires strictly increasing x_data.");
    }

    const double sig = hPrev / hTot;
    const double p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    const double dd =
        (ydata[i + 1] - ydata[i]) / hNext - (ydata[i] - ydata[i - 1]) / hPrev;
    u[i] = (6.0 * dd / hTot - sig * u[i - 1]) / p;
  }

  for (size_t k = n - 1; k-- > 0;) {
    y2[k] = y2[k] * y2[k + 1] + u[k];
  }
  return y2;
}

double CustomU::interpolateSpline(double z, const std::vector<double> &xdata,
                                  const std::vector<double> &ydata,
                                  const std::vector<double> &y2data) const {
  if (xdata.size() < 3 || ydata.size() != xdata.size() ||
      y2data.size() != xdata.size()) {
    throw std::invalid_argument(
        "Spline interpolation requires aligned x/y/y2 data with at least 3 "
        "points.");
  }

  if (z <= xdata.front()) {
    return ydata.front();
  }
  if (z >= xdata.back()) {
    return ydata.back();
  }

  auto upper = std::lower_bound(xdata.begin(), xdata.end(), z);
  size_t khi = static_cast<size_t>(std::distance(xdata.begin(), upper));
  size_t klo = khi - 1;
  const double h = xdata[khi] - xdata[klo];
  if (h <= 0.0) {
    throw std::invalid_argument(
        "Spline interpolation requires strictly increasing x_data.");
  }

  const double wa = (xdata[khi] - z) / h;
  const double wb = (z - xdata[klo]) / h;
  return wa * ydata[klo] + wb * ydata[khi] +
         ((wa * wa * wa - wa) * y2data[klo] +
          (wb * wb * wb - wb) * y2data[khi]) *
             (h * h) / 6.0;
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
  return diff_data;
}

void OSSolver::mapToStandardRegion() {
  // Map the flow profile to the standard region [-1, 1]
  CustomU *profile = dynamic_cast<CustomU *>(Uprof);
  if (!profile) {
    throw std::runtime_error("Uprofile is not of type CustomU");
  }

  a = profile->x_data[0];
  b = profile->x_data[profile->x_data.size() - 1];

  for (uint i = 0; i < profile->x_data.size(); i++) {
    profile->x_data[i] = profile->mapToStandardRegion(profile->x_data[i]);
  }
  profile->refreshInterpolationData();
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

void OSSolver::setArrays() {
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
      dN[i][j] = getdN(i + 1, z) * jacobian[j];
      d2N[i][j] = getd2N(i + 1, z) * jacobian[j] * jacobian[j] +
                  getdN(i + 1, z) * djacobian[j];
      if (i == 0) {
        U[j] = Uprof->getU(z);
        d2U[j] = Uprof->getd2U(z);
      }
    }
  }
}

void OSSolver::buildIntegralBlocks() {
  RealMatrix Nmat(dimVS, numQuadPoints);
  RealMatrix dNmat(dimVS, numQuadPoints);
  RealMatrix d2Nmat(dimVS, numQuadPoints);
  RealVector weights(numQuadPoints);
  RealVector weightsU(numQuadPoints);
  RealVector weightsd2U(numQuadPoints);

  for (uint i = 0; i < dimVS; ++i) {
    for (uint k = 0; k < numQuadPoints; ++k) {
      Nmat(i, k) = N[i][k];
      dNmat(i, k) = dN[i][k];
      d2Nmat(i, k) = d2N[i][k];
    }
  }

  for (uint k = 0; k < numQuadPoints; ++k) {
    const double w = gaussWeights[k] / jacobian[k];
    weights(k) = w;
    weightsU(k) = w * U[k];
    weightsd2U(k) = w * d2U[k];
  }

  const auto weightedProduct = [](const RealMatrix &lhs, const RealVector &w,
                                  const RealMatrix &rhs) -> RealMatrix {
    return (lhs.array().rowwise() * w.transpose().array()).matrix() *
           rhs.transpose();
  };

  T1 = weightedProduct(d2Nmat, weights, d2Nmat).cast<complex>();
  T2 = weightedProduct(d2Nmat, weights, Nmat).cast<complex>();
  T3 = weightedProduct(Nmat, weights, Nmat).cast<complex>();
  T4 = weightedProduct(Nmat, weightsd2U, Nmat).cast<complex>();
  T5 = weightedProduct(d2Nmat, weightsU, Nmat).cast<complex>();
  T6 = weightedProduct(Nmat, weightsU, Nmat).cast<complex>();
  T7 = weightedProduct(dNmat, weights, Nmat).cast<complex>();
  T8 = weightedProduct(d2Nmat, weightsU, dNmat).cast<complex>();
  T9 = weightedProduct(dNmat, weightsU, Nmat).cast<complex>();
}

// Build the A and B matrices for the generalized eigenvalue problem
void OSSolver::buildMatricesTemporal() {
  const complex I(0.0, 1.0);
  const complex alpha = var;

  A = T1 - 2.0 * k2 * T2 + k2 * k2 * T3 + I * alpha * re * T4 -
      I * alpha * re * T5 + I * alpha * k2 * re * T6;
  B = -I * re * (T2 - k2 * T3);
}

// Build the A and B matrices for the generalized eigenvalue problem
void OSSolver::buildMatricesSpatial() {
  const complex I(0.0, 1.0);
  const complex beta2 = beta * beta;
  const complex omega = var;
  const Eigen::Index dim = static_cast<Eigen::Index>(dimVS);

  const Matrix R0m = 1.0 / re * T1 - (2.0 / re * beta2 - I * omega) * T2 -
                     (I * omega * beta2 - 1.0 / re * beta2 * beta2) * T3;
  const Matrix R1m = I * T4 - I * T5 + I * beta2 * T6 -
                     (2.0 * I * omega - 4.0 * beta2 / re) * T7 + 4.0 / re * T8;
  const Matrix R2m = 4.0 / re * T2 + 2.0 * I * T9;

  A = Matrix::Zero(2 * dim, 2 * dim);
  B = Matrix::Zero(2 * dim, 2 * dim);

  A.topLeftCorner(dim, dim) = -R1m;
  A.topRightCorner(dim, dim) = -R0m;
  B.topLeftCorner(dim, dim) = R2m;
  A.bottomLeftCorner(dim, dim) = Matrix::Identity(dim, dim);
  B.bottomRightCorner(dim, dim) = Matrix::Identity(dim, dim);
}

Eigen::VectorXcd
OSSolver::computeEVec(const Eigen::VectorXcd &eigenvector_coeffs,
                      const std::string branch) const {
  const Eigen::Index dim = static_cast<Eigen::Index>(dimVS);
  const Eigen::Index expectedSize = (branch == BRANCH_SPATIAL) ? 2 * dim : dim;
  if (eigenvector_coeffs.size() < expectedSize) {
    throw std::invalid_argument("Eigenvector coefficient size is inconsistent "
                                "with the current branch.");
  }
  const Eigen::Index offset = (branch == BRANCH_SPATIAL) ? dim : 0;

  Eigen::VectorXcd eigenvector(numQuadPoints);
  for (uint j = 0; j < numQuadPoints; j++) {
    eigenvector[j] = 0.0;
    for (uint i = 0; i < dimVS; i++) {
      const Eigen::Index coeffIndex = static_cast<Eigen::Index>(i) + offset;
      eigenvector[j] += eigenvector_coeffs[coeffIndex] * N[i][j];
    }
  }

  return eigenvector;
}

// Solve Ax = lambda Bx using QZ (LAPACK zggev), avoiding explicit B^{-1}A.
EigenSolution OSSolver::solve() const {
  if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
    throw std::runtime_error("Generalized eigenvalue matrices must be square "
                             "and have matching dimensions.");
  }

  const Eigen::Index nEigen = A.rows();
  if (nEigen <= 0) {
    throw std::runtime_error("Generalized eigenvalue problem has zero size.");
  }
  if (nEigen > static_cast<Eigen::Index>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("Generalized eigenvalue problem is too large for "
                             "LAPACK integer interface.");
  }

  Matrix Awork = A;
  Matrix Bwork = B;
  Vector alphaVec(nEigen);
  Vector betaVec(nEigen);
  Matrix VR =
      Matrix::Zero(nEigen, nEigen); // Right eigenvectors (columns of VR)

  char jobvl = 'N'; // No left eigenvectors needed
  char jobvr = 'V'; // Compute right eigenvectors
  int n = static_cast<int>(nEigen);
  int lda = n;    // Leading dimension of A
  int ldb = n;    // Leading dimension of B
  int ldvl = 1;   // Leading dimension of VL (not used since jobvl = 'N')
  int ldvr = n;   // Leading dimension of VR
  int info = 0;   // Output info from zggev
  int lwork = -1; // Workspace query, means "don’t solve yet—tell me how much
                  // workspace memory you need"
  std::complex<double> workQuery = 0.0;
  std::vector<double> rwork(static_cast<size_t>(8 * n),
                            0.0); // Real workspace for zggev

  zggev_(&jobvl, &jobvr, &n, Awork.data(), &lda, Bwork.data(), &ldb,
         alphaVec.data(), betaVec.data(), nullptr, &ldvl, VR.data(), &ldvr,
         &workQuery, &lwork, rwork.data(), &info);
  if (info != 0) {
    throw std::runtime_error("LAPACK workspace query for generalized "
                             "eigensolver failed.");
  }

  lwork = std::max(1, static_cast<int>(std::real(workQuery)));
  std::vector<std::complex<double>> work(static_cast<size_t>(lwork), 0.0);
  Awork = A;
  Bwork = B;

  zggev_(&jobvl, &jobvr, &n, Awork.data(), &lda, Bwork.data(), &ldb,
         alphaVec.data(), betaVec.data(), nullptr, &ldvl, VR.data(), &ldvr,
         work.data(), &lwork, rwork.data(), &info);
  if (info != 0) {
    throw std::runtime_error("LAPACK generalized eigensolver failed.");
  }

  Vector eigenvalues(nEigen);
  const double betaTol = 1e-12;
  for (Eigen::Index i = 0; i < nEigen; ++i) {
    if (std::abs(betaVec[i]) <= betaTol) {
      const double inf = std::numeric_limits<double>::infinity();
      eigenvalues[i] = complex(inf, inf);
    } else {
      eigenvalues[i] = alphaVec[i] / betaVec[i];
    }
  }

  return EigenSolution{eigenvalues, VR};
}
