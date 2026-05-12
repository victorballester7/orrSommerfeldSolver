#ifndef PP_HPP
#define PP_HPP

#include "config.hpp"
#include "solver.hpp"
#include <complex>
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>

using complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

class PostProcess {
  inline static const std::string eValsFilename = "data/eigenvalues.dat";
  inline static const std::string eVecsFilename = "data/eigenvector.dat";

  std::vector<complex> eigenvalues;
  Matrix eigenvectors_coeffs;
  const size_t num_ev = 20;
  std::vector<size_t> evals_idx;
  size_t idx_target;
  complex rescaleFactor = 1;
  std::string evLabel;

  std::complex<double> rescaleFactorEV(const std::string branch,
                                       const complex var,
                                       const bool use_c) const;
  std::vector<size_t> getRankedFiniteEVIndices(const std::string branch,
                                               const std::string problem) const;

  size_t getIdxMostUnstableEV(bool useTargetEV,
                              std::complex<double> targetEV) const;

public:
  PostProcess(const EigenSolution &_eig, const Config &config) {
    eigenvalues =
        std::vector<complex>(_eig.eigenvalues.data(),
                             _eig.eigenvalues.data() + _eig.eigenvalues.size());
    // rescale eigenvalues to use c if possible
    rescaleFactor = rescaleFactorEV(config.branch, config.var, config.use_c);

    if (config.branch == BRANCH_TEMPORAL && std::abs(config.var) > 1e-10) {
      // rescale EVs by alpha to get c (for better post-processing)
      for (auto &ev : eigenvalues) {
        ev *= 1. / config.var;
      }
    }

    eigenvectors_coeffs = _eig.eigenvectors;

    evals_idx = getRankedFiniteEVIndices(config.branch, config.problem);

    for (auto &ev : eigenvalues) {
      if (config.branch == BRANCH_TEMPORAL && std::abs(config.var) > 1e-10) {
        ev *= config.var; // rescale back to omega
      }
      ev *= rescaleFactor;
    }
    evLabel = Config::getEVlabel(config.branch, config.use_c);
    idx_target = getIdxMostUnstableEV(config.useTargetEV, config.targetEV);
  }

  complex getTargetEVal() const;
  Vector getTargetEVec() const;

  // Print the spectrum
  void printSpectrum() const;

  // Write the spectrum to a file
  void writeToFile(const std::vector<double> &y, const Vector &eigenvector,
                   const Uprofile *Uprof) const;
  static void writeToFile(const std::vector<complex> &vars,
                          const std::vector<complex> &evals, Config &config) {
    // eigenvalues
    std::ofstream file(eValsFilename);
    if (!file.is_open()) {
      std::cerr << "Error opening file: " << eValsFilename << std::endl;
      return;
    }

    std::string evLabel = Config::getEVlabel(config.branch, config.use_c);
    std::string varLabel = Config::getVarlabel(config.branch);
    file << "# Re(" << evLabel << ")   Im(" << evLabel << ")   Re(" << varLabel
         << ")   Im(" << varLabel << ")" << std::endl;
    for (uint i = 0; i < evals.size(); i++) {
      file << evals[i].real() << " " << evals[i].imag() << " " << vars[i].real()
           << " " << vars[i].imag() << std::endl;
    }

    file.close();
  };
};
#endif // PP_HPP
