#include "../include/pp.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>

void PostProcess::rescaleEV(complex &ev, const std::string branch,
                            const complex var, const bool use_c) const {
  // Rescale eigenvalue by the factor alpha
  double tol_var = 1e-10;
  if (branch == BRANCH_TEMPORAL && abs(var) > tol_var && use_c) {
    // lambda is omega (so config.var = alpha), so we want c
    // return 1. / var;
    ev = ev / var;
  } else if (branch == BRANCH_SPATIAL && abs(var) > tol_var && use_c) {
    // lambda is alpha (so config.var = omega), so we want c
    ev = ev / var;
  }
  // else {
  //   // return 1.;
  // }
}

std::vector<size_t>
PostProcess::getRankedFiniteEVIndices(const std::string branch,
                                      const std::string problem) const {
  const size_t n = eigenvalues.size();
  std::vector<size_t> indices;
  indices.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    if (std::isfinite(eigenvalues[i].real()) &&
        std::isfinite(eigenvalues[i].imag())) {
      indices.push_back(i);
    }
  }

  if (branch == BRANCH_TEMPORAL) {

    if (problem == PB_POISEUILLE || problem == PB_COUETTE) {

      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        if (eigenvalues[lhs].imag() != eigenvalues[rhs].imag()) {
          return eigenvalues[lhs].imag() > eigenvalues[rhs].imag();
        }
        return eigenvalues[lhs].real() > eigenvalues[rhs].real();
      });

    } else if (problem == PB_BOUNDARY_LAYER || problem == PB_CUSTOM) {

      std::vector<size_t> filtered;
      filtered.reserve(indices.size());

      for (size_t idx : indices) {
        if (eigenvalues[idx].real() <= 0.8) {
          filtered.push_back(idx);
        }
      }

      indices.swap(filtered);

      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        if (eigenvalues[lhs].imag() != eigenvalues[rhs].imag()) {
          return eigenvalues[lhs].imag() > eigenvalues[rhs].imag();
        }
        return eigenvalues[lhs].real() < eigenvalues[rhs].real();
      });

    } else {
      throw std::invalid_argument("Unsupported temporal problem type: " +
                                  problem);
    }

  } else { // spatial branch
    if (problem == PB_POISEUILLE) {
      std::vector<size_t> filtered;
      filtered.reserve(indices.size());

      for (size_t idx : indices) {
        if (eigenvalues[idx].imag() <= -1 || eigenvalues[idx].imag() >= 2) {
          continue;
        }
        filtered.push_back(idx);
      }

      indices.swap(filtered);

      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        if (eigenvalues[lhs].imag() != eigenvalues[rhs].imag()) {
          return eigenvalues[lhs].imag() < eigenvalues[rhs].imag();
        }
        return eigenvalues[lhs].real() > eigenvalues[rhs].real();
      });
    } else if (problem == PB_COUETTE) {
      std::vector<size_t> filtered;
      filtered.reserve(indices.size());

      for (size_t idx : indices) {
        if (eigenvalues[idx].imag() < 0 || eigenvalues[idx].real() < 0) {
          continue;
        }
        filtered.push_back(idx);
      }

      indices.swap(filtered);

      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        if (eigenvalues[lhs].imag() != eigenvalues[rhs].imag()) {
          return eigenvalues[lhs].imag() < eigenvalues[rhs].imag();
        }
        return eigenvalues[lhs].real() < eigenvalues[rhs].real();
      });
    } else if (problem == PB_BOUNDARY_LAYER || problem == PB_CUSTOM) {
      std::vector<size_t> filtered;
      filtered.reserve(indices.size());

      for (size_t idx : indices) {
        if (eigenvalues[idx].imag() <= -1 || eigenvalues[idx].imag() >= 1 || eigenvalues[idx].real() <= 1.5 || eigenvalues[idx].real() >= 20) {
          continue;
        }
        filtered.push_back(idx);
      }

      indices.swap(filtered);

      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        if (eigenvalues[lhs].imag() != eigenvalues[rhs].imag()) {
          return eigenvalues[lhs].imag() < eigenvalues[rhs].imag();
        }
        return eigenvalues[lhs].real() > eigenvalues[rhs].real();
      });
    }
  }

  std::vector<size_t> first_indices;

  const size_t count = std::min(num_ev, indices.size());
  first_indices.assign(indices.begin(), indices.begin() + count);

  return first_indices;
}
size_t PostProcess::getIdxMostUnstableEV(bool useTargetEV,
                                         std::complex<double> targetEV) const {
  if (!useTargetEV) {
    return evals_idx[0];
  }

  size_t currentIdx = 0;
  double minDistance = std::abs(eigenvalues[currentIdx] - targetEV);

  for (size_t i = 1; i < eigenvalues.size(); i++) {
    double distance = std::abs(eigenvalues[i] - targetEV);
    if (distance < minDistance) {
      minDistance = distance;
      currentIdx = i;
    }
  }

  return currentIdx;
}

complex PostProcess::getTargetEVal() const { return eigenvalues[idx_target]; }

Eigen::VectorXcd PostProcess::getTargetEVec() const {
  return eigenvectors_coeffs.col(idx_target);
}

// Print the spectrum
void PostProcess::printSpectrum() const {
  for (auto idx : evals_idx) {
    std::cout << evLabel << " = " << eigenvalues[idx] << std::endl;
  }
}

void PostProcess::writeToFile(const std::vector<double> &y,
                              const Vector &eigenvector,
                              const Uprofile *Uprof) const {
  // eigenvalues
  std::ofstream file(eValsFilename);
  if (!file.is_open()) {
    throw std::runtime_error("Error opening file: " + eValsFilename);
  }
  file << "# Re(" << evLabel << ")   Im(" << evLabel << ")" << std::endl;
  for (const auto &lambda : eigenvalues) {
    file << lambda.real() << " " << lambda.imag() << std::endl;
  }
  file.close();

  // eigenvector
  std::ofstream file2(eVecsFilename);
  if (!file2.is_open()) {
    throw std::runtime_error("Error opening file: " + eVecsFilename);
  }

  file2 << "# y   Re(v)   Im(v)" << std::endl;
  std::vector<std::tuple<double, double, double>> rows;
  rows.reserve(y.size());
  for (size_t i = 0; i < y.size(); i++) {
    rows.emplace_back(Uprof->mapToPhysicalRegion(y[i]), eigenvector[i].real(),
                      eigenvector[i].imag());
  }
  std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
    return std::get<0>(a) < std::get<0>(b);
  });
  for (const auto &row : rows) {
    file2 << std::get<0>(row) << " " << std::get<1>(row) << " "
          << std::get<2>(row) << std::endl;
  }

  file2.close();
}
