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
  std::vector<complex> eigenvalues;
  Matrix eigenvectors;

  std::complex<double> rescaleEV(const std::complex<double> &lambda) const;
  std::complex<double>
  rescaleEVForSelection(const std::complex<double> &lambda) const;
  static bool isFiniteComplex(const complex &z);
  std::vector<size_t> getRankedFiniteEigenvalueIndices() const;

  Eigen::VectorXcd getMostUnstableEigenvector() const;
  // std::complex<double> blasiusScaling(const std::complex<double> &lambda)
  // const;

public:
  PostProcess(const EigenSolution &_eig) {
    eigenvalues =
        std::vector<complex>(_eig.eigenvalues.data(),
                             _eig.eigenvalues.data() + _eig.eigenvalues.size());
    eigenvectors = _eig.eigenvectors;
  }

  PostProcess(std::vector<complex> &_eigenvalues) : eigenvalues(_eigenvalues) {}
  complex getMostUnstableEigenvalue(bool useTargetEV, std::complex<double> targetEV) const;
  complex getMostUnstableEigenvalueNotScaled(bool useTargetEV, std::complex<double> targetEV) const;

  // Print the spectrum
  void printSpectrum(std::string evLabel) const;

  // Write the spectrum to a file
  void writeToFile(const OSSolver &solver) const;
  void writeToFile(const std::vector<complex> &vars) const;
};
#endif // PP_HPP
