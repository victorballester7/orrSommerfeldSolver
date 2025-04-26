#ifndef OSSOLVER_HPP
#define OSSOLVER_HPP

#include "read_conf.hpp" // Include the configuration header
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>

// Define complex number type
using complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

class Uprofile {
public:
  double jacobian = 1.0;

  // Flow profile: U(z) for plane Poiseuille flow
  virtual double getU(double z) = 0;

  // Flow profile second derivative: U''(z)
  virtual double getd2U(double z) = 0;
};

// Plane Poiseuille flow profile
class Poiseuille : public Uprofile {
public:
  Poiseuille() {}

  double getU(double z) override {
    return 1.0 - z * z; // Plane Poiseuille flow
  }

  double getd2U(__attribute__((unused)) double z) override { return -2.0; }
};

// Plane Couette flow profile
class Couette : public Uprofile {
public:
  Couette() {}

  double getU(double z) override {
    return z; // Plane Couette flow
  }

  double getd2U(__attribute__((unused)) double z) override { return 0.0; }
};

// Blasius flow profile
class CustomU : public Uprofile {
private:
  size_t len_data;

  void readFromFile(const std::string &filename, std::vector<double> &xdata,
                    uint colX, std::vector<double> &data, uint colY,
                    uint numSkipHeaderLines);

  double interpolate(double z, const std::vector<double> &xdata,
                     const std::vector<double> &ydata);

  std::vector<double> diffData(const std::vector<double> &xdata,
                               const std::vector<double> &ydata);

  double setJacobian() {
    // Compute the Jacobian of the mapping [a, b] -> [-1, 1]
    double a = x_data[0];
    double b = x_data.back();
    return 2.0 / (b - a);
  }

public:
  std::vector<double> x_data;
  std::vector<double> u_data;
  std::vector<double> du_data;
  std::vector<double> d2u_data;

  CustomU(std::string filename, uint colX, uint colY, uint numSkipHeaderLines) {
    readFromFile(filename, x_data, colX, u_data, colY, numSkipHeaderLines);
    len_data = u_data.size();
    du_data = diffData(x_data, u_data);
    d2u_data = diffData(x_data, du_data);
    jacobian = setJacobian();
  }

  double getU(double z) override { return interpolate(z, x_data, u_data); }

  double getd2U(double z) override { return interpolate(z, x_data, d2u_data); }
};

class OSSolver {
private:
  uint p;          // Polynomial degree
  uint dimVS;      // Dimension of the vector space of basis functions
  double re;       // Reynolds number
  complex var;   // Wavenumber
  complex beta;    // Wavenumber
  complex k2;      // Square of the wavenumber
  double jacobian; // Jacobian of the mapping [a, b] -> [-1, 1]

  uint numQuadPoints; // Number of quadrature points
  std::vector<double> gaussPoints;
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
  double getL(uint n, double x);

  // Internal shape function N_i(z) as defined in equation (10)
  double getN(uint i, double z);

  // Second derivative of shape function
  double getdN(uint i, double z);

  // Second derivative of shape function
  double getd2N(uint i, double z);

  // Determine the appropriate number of quadrature points
  uint setNumQuadPoints() {
    // The basis function have degree p. At most we are multiplying
    // U * phi * phi
    // 2q - 1 = p + p + degU
    // q = p + (degU + 1) / 2
    // For Poiseuille degU = 2
    return p + 2;
  }

  void setGaussPointsWeights();

  uint getDimVectorSpace() {
    // The dimension of the vector space of basis functions is p - 3
    return p - 3;
  }

  void mapToStandardRegion(CustomU &Uprofile);

  void setFunctions(Uprofile &Uprofile);

  // Matrices for temporal branch
  complex Los(uint i, uint j);
  complex M(uint i, uint j);

  // Matrices for spatial branch
  complex R0(uint i, uint j);
  complex R1(uint i, uint j);
  complex R2(uint i, uint j);

  // Build the A and B matrices for the generalized eigenvalue problem
  void buildMatricesTemporal();
  void buildMatricesSpatial();

public:
  OSSolver(Config &config)
      : p(config.p), re(config.re), var(config.var), beta(config.beta),
        k2(config.k2) {

    assert(p > 3); // Ensure polynomial degree is greater than 3

    // Determine appropriate number of quadrature points
    numQuadPoints = setNumQuadPoints();

    // Initialize Gauss points and weights
    setGaussPointsWeights();

    // Set the dimension of the vector space of basis functions
    dimVS = getDimVectorSpace();

    // Initialize shape functions and their second derivatives
    Uprofile *Uprofile;

    if (config.problem == PB_POISEUILLE) {
      Uprofile = new Poiseuille();
    } else if (config.problem == PB_COUETTE) {
      Uprofile = new Couette();
    } else if (config.problem == PB_BOUNDARY_LAYER) {
      uint colX = 1;
      uint colY = 2;
      uint numSkipHeaderLines = 3;
      Uprofile =
          new CustomU("data/blasius.dat", colX, colY, numSkipHeaderLines);
      mapToStandardRegion(*dynamic_cast<CustomU *>(Uprofile));
    } else if (config.problem == PB_CUSTOM) {
      Uprofile = new CustomU(config.filenameUprofile, config.colX, config.colY,
                             config.numSkipHeaderLines);
      mapToStandardRegion(*dynamic_cast<CustomU *>(Uprofile));
    } else {
      throw std::invalid_argument(
          "Invalid profile type or not enough arguments provided");
    }

    jacobian = Uprofile->jacobian;

    setFunctions(*Uprofile);
  }

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

  // Solve the generalized eigenvalue problem
  std::vector<complex> solve();
};

#endif // OSSOLVER_HPP
