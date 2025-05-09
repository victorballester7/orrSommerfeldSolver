#ifndef PP_HPP
#define PP_HPP

#include "config.hpp"
#include "solver.hpp"
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>

using complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

class PostProcess {
  Config &config;
  Eigen::ComplexEigenSolver<Matrix> *eig;
  std::vector<complex> eigenvalues;
  Matrix eigenvectors;

  std::complex<double> rescaleEV(const std::complex<double> &alpha,
                                 const std::complex<double> &lambda) const;

  Eigen::VectorXcd getMostUnstableEigenvector() const;

public:
  PostProcess(Config &_config, Eigen::ComplexEigenSolver<Matrix> &_eig)
      : config(_config), eig(&_eig) {

    eigenvalues = std::vector<complex>(eig->eigenvalues().data(),
                                       eig->eigenvalues().data() +
                                           eig->eigenvalues().size());
    eigenvectors = eig->eigenvectors();
  }

  PostProcess(Config &_config, std::vector<complex> &_eigenvalues)
      : config(_config), eig(nullptr), eigenvalues(_eigenvalues) {}
  complex getMostUnstableEigenvalue() const;
  complex getMostUnstableEigenvalueNotScaled() const;

  // Print the spectrum
  void printSpectrum() const;

  // Write the spectrum to a file
  void writeToFile(const OSSolver &solver) const;
  void writeToFile(const std::vector<complex> &vars) const;
  // Plot the spectrum
  // void plotSpectrum();
};
#endif // PP_HPP
