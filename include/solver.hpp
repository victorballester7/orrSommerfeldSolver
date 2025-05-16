#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "config.hpp" // Include the configuration header
#include "flowProfiles.hpp"
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>

// Define complex number type
using complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

class OSSolver {
private:
  uint p;          // Polynomial degree
  uint dimVS;      // Dimension of the vector space of basis functions
  double re;       // Reynolds number
  complex var;     // Wavenumber
  complex beta;    // Wavenumber
  complex k2;      // Square of the wavenumber
  double a = -1.0; // Left boundary of the physical region
  double b = 1.0;  // Right boundary of the physical region

  std::vector<double> jacobian;  // Jacobian of the mapping [a, b] -> [-1, 1]
  std::vector<double> djacobian; // Derivative of the Jacobian

  std::vector<double> gaussWeights;

  // matrix of shape functions (row ith contains N_i(z_j) for every j)
  std::vector<std::vector<double>> N;   // Shape functions
  std::vector<std::vector<double>> dN;  // First derivatives of shape functions
  std::vector<std::vector<double>> d2N; // Second derivatives of shape functions

  std::vector<double> U;   // Flow profile
  std::vector<double> d2U; // Second derivatives of flow profile

  // Matrices for the eigenvalue problem Ax = λBx
  Matrix A;
  Matrix B;

  // Evaluate Legendre polynomial of degree n at x using recursion formula
  double getL(uint n, double x) const;

  // Internal shape function N_i(z) as defined in equation (10)
  double getN(uint i, double z) const;

  // Second derivative of shape function
  double getdN(uint i, double z) const;

  // Second derivative of shape function
  double getd2N(uint i, double z) const;

  // Determine the appropriate number of quadrature points
  uint getNumQuadPoints() const {
    // The basis function have degree p. At most we are multiplying
    // U * phi * phi
    // 2q - 1 = p + p + degU
    // q = p + (degU + 1) / 2
    // For Poiseuille degU = 2
    return p + 2;
  }

  void setGaussPointsWeights();

  uint getDimVectorSpace() const {
    // The dimension of the vector space of basis functions is p - 3
    return p - 3;
  }

  void mapToStandardRegion();

  void setFunctions();

  // Matrices for temporal branch
  complex Los(uint i, uint j) const;
  complex M(uint i, uint j) const;

  // Matrices for spatial branch
  complex R0(uint i, uint j) const;
  complex R1(uint i, uint j) const;
  complex R2(uint i, uint j) const;

  // Build the A and B matrices for the generalized eigenvalue problem
  void buildMatricesTemporal();
  void buildMatricesSpatial();

public:
  uint numQuadPoints; // Number of quadrature points
  std::vector<double> gaussPoints;
  Uprofile *Uprof;

  OSSolver(Config &config)
      : p(config.p), re(config.re), var(config.var), beta(config.beta),
        k2(config.k2) {

    assert(p > 3); // Ensure polynomial degree is greater than 3

    // Determine appropriate number of quadrature points
    numQuadPoints = getNumQuadPoints();

    // Initialize Gauss points and weights
    setGaussPointsWeights();

    // Set the dimension of the vector space of basis functions
    dimVS = getDimVectorSpace();

    // Initialize shape functions and their second derivatives

    if (config.problem == PB_POISEUILLE) {
      Uprof = new Poiseuille();
    } else if (config.problem == PB_COUETTE) {
      Uprof = new Couette();
    } else if (config.problem == PB_BOUNDARY_LAYER) {
      uint colX = 1;
      uint colY = 2;
      uint numSkipHeaderLines = 3;
      Uprof =
          new CustomU("data/blasius.dat", colX, colY, numSkipHeaderLines);
    } else if (config.problem == PB_CUSTOM) {
      Uprof = new CustomU(config.filenameUprofile, config.colX, config.colY,
                             config.numSkipHeaderLines);
    } else {
      throw std::invalid_argument(
          "Invalid profile type or not enough arguments provided");
    }

    jacobian = Uprof->getJacobian(gaussPoints);
    djacobian = Uprof->getDiffJacobian(gaussPoints);

    // if Uprofile is a CustomU, map to standard region
    if (config.problem == PB_BOUNDARY_LAYER || config.problem == PB_CUSTOM) {
      mapToStandardRegion();
    }

    setFunctions();
  }

  double getYPhysicalRegion(const OSSolver &solver, uint i) const;

  void buildMatrices(std::string branch) {
    if (branch == BRANCH_TEMPORAL) {
      buildMatricesTemporal();
    } else if (branch == BRANCH_SPATIAL) {
      buildMatricesSpatial();
    }
  }

  void setVar(complex var_, std::string branch) {
    var = var_;
    if (branch == BRANCH_TEMPORAL) {
      k2 = var * var + beta * beta;
    }
  }

  Eigen::VectorXcd
  computeEigenvector(const Eigen::VectorXcd &eigenvector_coeffs) const;

  // Solve the generalized eigenvalue problem
  Eigen::ComplexEigenSolver<Matrix> solve() const;
};

#endif // SOLVER_HPP
