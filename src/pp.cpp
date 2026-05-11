#include "../include/pp.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace {
double quantile(std::vector<double> values, double q) {
  if (values.empty()) {
    return 0.0;
  }
  q = std::clamp(q, 0.0, 1.0);
  std::sort(values.begin(), values.end());

  const double pos = q * static_cast<double>(values.size() - 1);
  const size_t lo = static_cast<size_t>(std::floor(pos));
  const size_t hi = static_cast<size_t>(std::ceil(pos));
  if (lo == hi) {
    return values[lo];
  }
  const double t = pos - static_cast<double>(lo);
  return values[lo] + (values[hi] - values[lo]) * t;
}

double median(const std::vector<double> &values) {
  return quantile(values, 0.5);
}

double robustScale(const std::vector<double> &values) {
  if (values.empty()) {
    return 1.0;
  }

  const double med = median(values);
  std::vector<double> absDev;
  absDev.reserve(values.size());
  for (const double value : values) {
    absDev.push_back(std::abs(value - med));
  }

  double scale = 1.4826 * median(absDev); // MAD to std-dev scale
  if (scale < 1e-12) {
    const double iqr = quantile(values, 0.75) - quantile(values, 0.25);
    scale = iqr / 1.349;
  }
  return std::max(scale, 1e-10);
}

double adaptiveLowerBound(const std::vector<double> &values,
                          double sigmaFactor = 3.5) {
  return median(values) - sigmaFactor * robustScale(values);
}

double adaptiveUpperBound(const std::vector<double> &values,
                          double sigmaFactor = 3.5) {
  return median(values) + sigmaFactor * robustScale(values);
}

double inferSpatialImagUpper(const std::vector<double> &absImagValues) {
  if (absImagValues.empty()) {
    return std::numeric_limits<double>::infinity();
  }
  if (absImagValues.size() < 8) {
    return adaptiveUpperBound(absImagValues, 2.0);
  }

  std::vector<double> logAbs;
  logAbs.reserve(absImagValues.size());
  for (const double value : absImagValues) {
    logAbs.push_back(std::log10(std::max(value, 1e-12)));
  }

  double c1 = quantile(logAbs, 0.25);
  double c2 = quantile(logAbs, 0.75);

  for (int iter = 0; iter < 25; ++iter) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    size_t n1 = 0;
    size_t n2 = 0;
    for (const double value : logAbs) {
      if (std::abs(value - c1) <= std::abs(value - c2)) {
        sum1 += value;
        ++n1;
      } else {
        sum2 += value;
        ++n2;
      }
    }
    if (n1 == 0 || n2 == 0) {
      break;
    }
    const double nextC1 = sum1 / static_cast<double>(n1);
    const double nextC2 = sum2 / static_cast<double>(n2);
    if (std::abs(nextC1 - c1) + std::abs(nextC2 - c2) < 1e-8) {
      c1 = nextC1;
      c2 = nextC2;
      break;
    }
    c1 = nextC1;
    c2 = nextC2;
  }

  if (c1 > c2) {
    std::swap(c1, c2);
  }

  const double fallback = adaptiveUpperBound(absImagValues, 2.0);
  if ((c2 - c1) < 0.4) {
    return fallback;
  }
  const double clusterSplit = std::pow(10.0, 0.5 * (c1 + c2));
  return std::min(clusterSplit, fallback);
}

double inferSpatialColumnHalfWidth(const std::vector<double> &absOffsets) {
  if (absOffsets.size() < 12) {
    return 0.0;
  }

  std::vector<double> positiveOffsets;
  positiveOffsets.reserve(absOffsets.size());
  for (const double value : absOffsets) {
    if (value > 1e-10) {
      positiveOffsets.push_back(value);
    }
  }
  if (positiveOffsets.size() < 12) {
    return 0.0;
  }

  std::vector<double> logOffsets;
  logOffsets.reserve(positiveOffsets.size());
  for (const double value : positiveOffsets) {
    logOffsets.push_back(std::log10(std::max(value, 1e-14)));
  }

  double c1 = quantile(logOffsets, 0.15);
  double c2 = quantile(logOffsets, 0.85);

  for (int iter = 0; iter < 25; ++iter) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    size_t n1 = 0;
    size_t n2 = 0;
    for (const double value : logOffsets) {
      if (std::abs(value - c1) <= std::abs(value - c2)) {
        sum1 += value;
        ++n1;
      } else {
        sum2 += value;
        ++n2;
      }
    }
    if (n1 == 0 || n2 == 0) {
      break;
    }
    const double nextC1 = sum1 / static_cast<double>(n1);
    const double nextC2 = sum2 / static_cast<double>(n2);
    if (std::abs(nextC1 - c1) + std::abs(nextC2 - c2) < 1e-8) {
      c1 = nextC1;
      c2 = nextC2;
      break;
    }
    c1 = nextC1;
    c2 = nextC2;
  }

  if (c1 > c2) {
    std::swap(c1, c2);
  }

  if ((c2 - c1) < 0.6) {
    return 0.0;
  }

  const double splitLog = 0.5 * (c1 + c2);
  const double split = std::pow(10.0, splitLog);
  size_t lowCount = 0;
  for (const double value : positiveOffsets) {
    if (value <= split) {
      ++lowCount;
    }
  }
  const size_t minLowCount = std::max<size_t>(
      10,
      static_cast<size_t>(0.25 * static_cast<double>(positiveOffsets.size())));
  if (lowCount < minLowCount) {
    return 0.0;
  }

  return split;
}
} // namespace

// std::complex<double>
// PostProcess::blasiusScaling(const std::complex<double> &lambda) const {
//   if (config.problem == PB_BOUNDARY_LAYER) {
//     // rescale quantities for Blasius flow because delta* is not 1, it is
//     // DELTASTAR_BLASIUS
//     return lambda * DELTASTAR_BLASIUS;
//   } else {
//     return lambda;
//   }
// }

std::complex<double>
PostProcess::rescaleEV(const std::complex<double> &lambda) const {
  // Rescale eigenvalue by the factor alpha
  double tol_alpha = 1e-10;
  if (config.branch == BRANCH_TEMPORAL && abs(config.var) > tol_alpha &&
      config.use_c) {
    return lambda /
           config.var; // lambda is omega (so config.var = alpha), so we want c
  } else {
    return lambda;
  }
}


bool PostProcess::isFiniteComplex(const complex &z) {
  return std::isfinite(z.real()) && std::isfinite(z.imag());
}

std::vector<size_t> PostProcess::getRankedFiniteEigenvalueIndices() const {
  std::vector<size_t> indices;
  indices.reserve(eigenvalues.size());

  std::vector<complex> scaled(eigenvalues.size(), complex(0.0, 0.0));
  for (size_t i = 0; i < eigenvalues.size(); ++i) {
    if (!isFiniteComplex(eigenvalues[i])) {
      continue;
    }
    scaled[i] = rescaleEVForSelection(eigenvalues[i]);
    if (isFiniteComplex(scaled[i])) {
      indices.push_back(i);
    }
  }

  if (indices.empty()) {
    throw std::runtime_error("No finite eigenvalues available.");
  }

  auto applyFilterWithFallback = [&](const auto &predicate) {
    std::vector<size_t> filtered;
    filtered.reserve(indices.size());
    for (const size_t idx : indices) {
      if (predicate(idx)) {
        filtered.push_back(idx);
      }
    }
    if (!filtered.empty()) {
      indices.swap(filtered);
    }
  };

  const bool isTemporal = (config.branch == BRANCH_TEMPORAL);
  const bool isBoundaryLike =
      (config.problem == PB_BOUNDARY_LAYER || config.problem == PB_CUSTOM);
  const bool isPoiseuille = (config.problem == PB_POISEUILLE);

  if (isTemporal && isBoundaryLike) {
    std::vector<double> ci;
    std::vector<double> cr;
    ci.reserve(indices.size());
    cr.reserve(indices.size());
    for (const size_t idx : indices) {
      ci.push_back(scaled[idx].imag());
      cr.push_back(scaled[idx].real());
    }

    const double ciLower = adaptiveLowerBound(ci, 3.0);
    const double crUpper = adaptiveUpperBound(cr, 3.0);
    applyFilterWithFallback([&](size_t idx) {
      return scaled[idx].imag() >= ciLower && scaled[idx].real() <= crUpper;
    });

    // Remove the dense continuous branch around c_r ~= 1 by splitting c_r into
    // two clusters and dropping the dominant upper-real cluster when present.
    if (indices.size() >= 20) {
      std::vector<double> crValues;
      crValues.reserve(indices.size());
      for (const size_t idx : indices) {
        crValues.push_back(scaled[idx].real());
      }

      double cLow = quantile(crValues, 0.20);
      double cHigh = quantile(crValues, 0.80);
      for (int iter = 0; iter < 30; ++iter) {
        double sumLow = 0.0;
        double sumHigh = 0.0;
        size_t nLow = 0;
        size_t nHigh = 0;
        for (const double crVal : crValues) {
          if (std::abs(crVal - cLow) <= std::abs(crVal - cHigh)) {
            sumLow += crVal;
            ++nLow;
          } else {
            sumHigh += crVal;
            ++nHigh;
          }
        }
        if (nLow == 0 || nHigh == 0) {
          break;
        }
        const double nextLow = sumLow / static_cast<double>(nLow);
        const double nextHigh = sumHigh / static_cast<double>(nHigh);
        if (std::abs(nextLow - cLow) + std::abs(nextHigh - cHigh) < 1e-10) {
          cLow = nextLow;
          cHigh = nextHigh;
          break;
        }
        cLow = nextLow;
        cHigh = nextHigh;
      }
      if (cLow > cHigh) {
        std::swap(cLow, cHigh);
      }

      const double split = 0.5 * (cLow + cHigh);
      const double separation = cHigh - cLow;
      size_t highCount = 0;
      for (const double crVal : crValues) {
        if (crVal >= split) {
          ++highCount;
        }
      }
      const size_t minHighCount = std::max<size_t>(
          12, static_cast<size_t>(0.35 * static_cast<double>(indices.size())));
      if (cHigh > 0.8 && separation > 0.15 && highCount >= minHighCount) {
        applyFilterWithFallback(
            [&](size_t idx) { return scaled[idx].real() < split; });
      }
    }
  } else if (!isTemporal) {
    std::vector<double> absAlphaI;
    absAlphaI.reserve(indices.size());
    for (const size_t idx : indices) {
      if (scaled[idx].real() > 0.0) {
        absAlphaI.push_back(std::abs(scaled[idx].imag()));
      }
    }
    if (absAlphaI.empty()) {
      for (const size_t idx : indices) {
        absAlphaI.push_back(std::abs(scaled[idx].imag()));
      }
    }

    double alphaIUpper = inferSpatialImagUpper(absAlphaI);
    if (isPoiseuille) {
      const double lowTail = quantile(absAlphaI, 0.05);
      const double cap = 1.2 * quantile(absAlphaI, 0.10);
      alphaIUpper = std::min(alphaIUpper, std::max(lowTail, cap));
    }
    if (isPoiseuille) {
      applyFilterWithFallback([&](size_t idx) {
        return scaled[idx].real() > 0.0 &&
               std::abs(scaled[idx].imag()) <= alphaIUpper;
      });
    } else {
      std::vector<double> alphaR;
      alphaR.reserve(indices.size());
      for (const size_t idx : indices) {
        alphaR.push_back(scaled[idx].real());
      }
      const double alphaRUpperRough = adaptiveUpperBound(alphaR, 3.0);
      applyFilterWithFallback([&](size_t idx) {
        return scaled[idx].real() > 0.0 &&
               scaled[idx].real() <= alphaRUpperRough &&
               std::abs(scaled[idx].imag()) <= alphaIUpper;
      });
    }

    if (isPoiseuille) {
      // Keep the most distinct spatial branch: largest alpha_r among the
      // low-|alpha_i| candidates.
      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        const complex &a = scaled[lhs];
        const complex &b = scaled[rhs];
        if (a.real() != b.real()) {
          return a.real() > b.real();
        }
        if (std::abs(a.imag()) != std::abs(b.imag())) {
          return std::abs(a.imag()) < std::abs(b.imag());
        }
        return a.imag() < b.imag();
      });
      return indices;
    } else {
      if (config.problem != PB_BOUNDARY_LAYER) {
        std::vector<double> alphaRCore;
        alphaRCore.reserve(indices.size());
        for (const size_t idx : indices) {
          alphaRCore.push_back(scaled[idx].real());
        }
        const double alphaRUpperRefined = adaptiveUpperBound(alphaRCore, 2.5);
        applyFilterWithFallback([&](size_t idx) {
          return scaled[idx].real() <= alphaRUpperRefined;
        });
      }

      if (config.problem == PB_BOUNDARY_LAYER && indices.size() >= 12) {
        const double columnCenter =
            config.var.real(); // alpha_r ~= omega_r column
        std::vector<double> absOffsets;
        absOffsets.reserve(indices.size());
        for (const size_t idx : indices) {
          absOffsets.push_back(std::abs(scaled[idx].real() - columnCenter));
        }

        const double columnHalfWidth = inferSpatialColumnHalfWidth(absOffsets);
        if (columnHalfWidth > 0.0) {
          const double rejectionHalfWidth = 2.0 * columnHalfWidth;
          applyFilterWithFallback([&](size_t idx) {
            return std::abs(scaled[idx].real() - columnCenter) >=
                   rejectionHalfWidth;
          });
        }
      }
    }
  }

  if (isTemporal) {
    std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
      const complex &a = scaled[lhs];
      const complex &b = scaled[rhs];
      if (a.imag() != b.imag()) {
        return a.imag() > b.imag();
      }
      if (isBoundaryLike && a.real() != b.real()) {
        return a.real() < b.real();
      }
      return a.real() > b.real();
    });
  } else {
    if (config.problem == PB_BOUNDARY_LAYER) {
      const double columnCenter = config.var.real();
      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        const complex &a = scaled[lhs];
        const complex &b = scaled[rhs];
        const double offA = std::abs(a.real() - columnCenter);
        const double offB = std::abs(b.real() - columnCenter);
        if (offA != offB) {
          return offA < offB; // closest mode off alpha_r = omega_r column
        }
        const double absAiA = std::abs(a.imag());
        const double absAiB = std::abs(b.imag());
        if (absAiA != absAiB) {
          return absAiA < absAiB;
        }
        if (a.real() != b.real()) {
          return a.real() > b.real();
        }
        return a.imag() < b.imag();
      });
    } else {
      std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs) {
        const complex &a = scaled[lhs];
        const complex &b = scaled[rhs];
        if (a.real() != b.real()) {
          return a.real() > b.real(); // BL/custom: prioritize alpha_r
        }
        if (a.imag() != b.imag()) {
          return a.imag() < b.imag(); // prefer larger -alpha_i
        }
        return std::abs(a.imag()) < std::abs(b.imag());
      });
    }
  }

  return indices;
}

size_t PostProcess::getIdxMostUnstableEVal(bool useTargetEV, std::complex<double> targetEV) const {
  if (eigenvalues.empty()) {
    throw std::runtime_error("No eigenvalues available.");
  }

  if (useTargetEV) {
    double minDistance = std::numeric_limits<double>::infinity();
    size_t bestIdx = 0;
    bool found = false;
    for (size_t idx = 0; idx < eigenvalues.size(); ++idx) {
      if (!isFiniteComplex(eigenvalues[idx])) { //necessary check due to solving of lambda = alpha/beta in zggev
        continue;
      }
      const complex lambdaCompare = rescaleEVForSelection(eigenvalues[idx]);
      // const complex lambdaCompare = eigenvalues[idx];
      if (!isFiniteComplex(lambdaCompare)) {
        continue;
      }
      const double distance = std::abs(lambdaCompare - targetEV);
      if (distance < minDistance) {
        minDistance = distance;
        bestIdx = idx;
        found = true;
      }
    }
    if (!found) {
      throw std::runtime_error(
          "No finite eigenvalues available for targetEV selection.");
    }
    return bestIdx;
  }

  const auto ranked = getRankedFiniteEigenvalueIndices();
  return ranked.front();
}

complex PostProcess::getMostUnstableEigenvalue(bool useTargetEV, std::complex<double> targetEV) const {
  complex max_eigenvalue = getMostUnstableEigenvalueNotScaled(useTargetEV, targetEV);

  return rescaleEV(max_eigenvalue);
}

Eigen::VectorXcd PostProcess::getMostUnstableEigenvector() const {
  // Find the eigenvalue with the largest imaginary part
  std::complex<double> max_eigenvalue = getMostUnstableEigenvalueNotScaled();

  // Find the index of that eigenvalue
  int index = -1;
  for (size_t i = 0; i < eigenvalues.size(); ++i) {
    if (std::abs(eigenvalues[i] - max_eigenvalue) < 1e-12) {
      index = static_cast<int>(i);
      break;
    }
  }

  if (index == -1) {
    throw std::runtime_error(
        "Most unstable eigenvalue not found in eigenvalue list.");
  }

  // Extract the corresponding eigenvector
  Eigen::VectorXcd vec = eigenvectors.col(index);

  return vec;
}

// Print the spectrum
void PostProcess::printSpectrum(const std::string evLabel) const {
  const size_t LIMIT = 20; // Limit output to first 20 eigenvalues
  const auto ranked = getRankedFiniteEigenvalueIndices();
  for (size_t i = 0; i < std::min(LIMIT, ranked.size()); i++) {
    const complex lambda = rescaleEV(eigenvalues[ranked[i]]);
    std::cout << evLabel << " = " << lambda.real() << " + " << lambda.imag()
              << "i" << std::endl;
  }
}

void PostProcess::writeEValsToFile(const std::string evLabel, const std::string filename) const {
std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  file << "# Re(" << evLabel << ")   Im(" << evLabel << ")" << std::endl;
  for (const auto &lambda : eigenvalues) {
    auto lambda_scaled = rescaleEV(lambda);
    file << lambda_scaled.real() << " " << lambda_scaled.imag() << std::endl;
  }
  file.close();
}

void PostProcess::writeEVecsToFile(const std::vector<double> &yPhysical, const std::string filename) const {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  file << "# y   Re(v)   Im(v)" << std::endl;
  Eigen::VectorXcd eigenvector = getMostUnstableEigenvector();
  std::vector<std::tuple<double, double, double>> rows;
  rows.reserve(eigenvector.size());
  for (int i = 0; i < eigenvector.size(); i++) {
    rows.emplace_back(yPhysical[i], eigenvector[i].real(), eigenvector[i].imag());
  }
  std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
    return std::get<0>(a) < std::get<0>(b);
  });
  for (const auto &row : rows) {
    file << std::get<0>(row) << " " << std::get<1>(row) << " "
         << std::get<2>(row) << std::endl;
  }
  file.close();
}


void PostProcess::writeToFile(const std::string evLabel,
                              const OSSolver &solver) const {
  std::string eValsFilename = "data/eigenvalues.dat";
  std::string eVecsFilename = "data/eigenvector.dat";

  // eigenvalues
  std::ofstream file(eValsFilename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << eValsFilename << std::endl;
    return;
  }
  file << "# Re(" << evLabel << ")   Im(" << evLabel << ")" << std::endl;
  for (const auto &lambda : eigenvalues) {
    auto lambda_scaled = rescaleEV(lambda);
    file << lambda_scaled.real() << " " << lambda_scaled.imag() << std::endl;
  }
  file.close();

  // eigenvector
  std::ofstream file2(eVecsFilename);
  if (!file2.is_open()) {
    std::cerr << "Error opening file: " << eVecsFilename << std::endl;
    return;
  }
  Eigen::VectorXcd eigenvector =
      solver.computeEigenvector(getMostUnstableEigenvector(), config.branch);

  file2 << "# y   Re(v)   Im(v)" << std::endl;
  std::vector<std::tuple<double, double, double>> rows;
  rows.reserve(solver.numQuadPoints);
  for (uint i = 0; i < solver.numQuadPoints; i++) {
    rows.emplace_back(solver.getYPhysicalRegion(solver, i),
                      eigenvector[i].real(), eigenvector[i].imag());
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

void PostProcess::writeToFile(const std::vector<complex> &var_complex) const {
  // eigenvalues
  std::ofstream file(config.filenameEigenvalues);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << config.filenameEigenvalues
              << std::endl;
    return;
  }

  std::string evLabel = config.getEVlabel();
  std::string varLabel = config.getVarlabel();
  file << "# Re(" << evLabel << ")   Im(" << evLabel << ")   Re(" << varLabel
       << ")   Im(" << varLabel << ")" << std::endl;
  for (uint i = 0; i < eigenvalues.size(); i++) {
    file << eigenvalues[i].real() << " " << eigenvalues[i].imag() << " "
         << var_complex[i].real() << " " << var_complex[i].imag() << std::endl;
  }

  file.close();
}

// void PostProcess::plotSpectrum() {
//   using namespace matplot;
//   // auto theta = linspace(0, 1, 500);
//   // auto x = transform(
//   //     theta, [&](double theta) { return exp(theta) * sin(100 * theta);
//   });
//   // auto y = transform(
//   //     theta, [&](double theta) { return exp(theta) * cos(100 * theta);
//   });

//   // auto s = scatter(x, y);
//   // s->marker_color("b");
//   // s->marker_face_color({0, .5, .5});
//   // At the beginning of your function
//   // show();
//   // Extract real and imaginary parts
//   std::vector<double> re, im;
//   for (const auto &eig : eigenvalues) {
//     re.push_back(eig.real()); // Real part
//     im.push_back(eig.imag()); // Imaginary part
//   }

//   // Create a scatter plot
//   // for (int i = 0; i < re.size(); i += chunk_size_points) {
//   //   auto re_aux = std::vector<double>(re.begin() + i,
//   //                         re.begin() + std::min(i + chunk_size_points,
//   //                         (int)re.size()));
//   //   auto im_aux = std::vector<double>(im.begin() + i,
//   //                         im.begin() + std::min(i + chunk_size_points,
//   //                         (int)im.size()));
//   //   auto s = scatter(re_aux, im_aux);
//   // }

//   auto s = scatter(re, im); // Plot with red color and size 10
//   s->marker_color("b");
//   s->marker_face_color({0, .5, .5});

//   // add grid
//   grid(true);

//   // Set axis limits
//   xlim({config.plot_lims.xmin, config.plot_lims.xmax});
//   ylim({config.plot_lims.ymin, config.plot_lims.ymax});

//   show();
// }
