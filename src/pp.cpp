#include "../include/pp.hpp"
#include <ostream>
#include <string>

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

std::complex<double>
PostProcess::blasiusScaling(const std::complex<double> &lambda) const {
  if (config.problem == PB_BOUNDARY_LAYER) {
    // rescale quantities for Blasius flow because delta* is not 1, it is
    // DELTASTAR_BLASIUS
    return lambda * DELTASTAR_BLASIUS;
  } else {
    return lambda;
  }
}

std::complex<double>
PostProcess::rescaleEV(const std::complex<double> &alpha,
                       const std::complex<double> &lambda) const {
  // Rescale eigenvalue by the factor alpha
  double tol_alpha = 1e-10;
  if (config.branch == BRANCH_TEMPORAL && abs(config.var) > tol_alpha &&
      config.use_c) {
    return lambda / alpha; // lambda is omega, so we want c, which doesn't
                           // need DELTASTAR rescaling
  } else {
    return blasiusScaling(lambda);
  }
}

complex PostProcess::getMostUnstableEigenvalueNotScaled() const {
  // Find the eigenvalue with the largest imaginary part
  complex max_eigenvalue = eigenvalues[0];
  uint numEVexamine = 20;
  // looks like the eigenvalue solver positions the target eigenvalue up in
  // the list
  for (uint i = 0; i < numEVexamine; i++) {
    if (eigenvalues[i].imag() > max_eigenvalue.imag()) {
      max_eigenvalue = eigenvalues[i];
    }
  }
  if (config.branch == BRANCH_TEMPORAL) {
    if (config.problem == PB_BOUNDARY_LAYER || config.problem == PB_CUSTOM) {
      double max_real = -1000;
      double min_real = 1000;

      for (const auto &lambda : eigenvalues) {
        if (lambda.real() > max_real) {
          max_real = lambda.real();
        }
        if (lambda.real() < min_real) {
          min_real = lambda.real();
        }
      }
      // printf("max_real = %f\n", max_real);
      // printf("min_real = %f\n", min_real);
      double tol_max_real = max_real - (max_real - min_real) * 0.3;

      // printf("tol_max_real = %f\n", tol_max_real);
      for (const auto &lambda : eigenvalues) {
        // std::cout << "eigenvalues[i].real() = " << lambda.real()
        //           << " + " << lambda.imag() << "i" << std::endl;
        // std::cout << "ev.imag > max_eigenvalue.imag() = "
        //           << (lambda.imag() > max_eigenvalue.imag())
        //           << std::endl;
        // std::cout << "ev.real < tol_max_real = "
        //           << (lambda.real() < tol_max_real) << std::endl;
        if ((lambda.imag() > max_eigenvalue.imag() ||
             max_eigenvalue.real() > tol_max_real ||
             max_eigenvalue.real() < 0) &&
            lambda.real() < tol_max_real && lambda.real() > 0) {
          // std::cout << "insideeee eigenvalues[i].real() = "
          //           << eigenvalues[i].real() << " + " <<
          //           eigenvalues[i].imag()
          //           << "i" << std::endl;
          max_eigenvalue = lambda;
        }
        //   std::cout << "max_eigenvalue = " << max_eigenvalue.real() << " +
        //   "
        //             << max_eigenvalue.imag() << "i" << std::endl;
      }
    } else {
      for (uint i = 0; i < numEVexamine; i++) {
        if (eigenvalues[i].imag() > max_eigenvalue.imag()) {
          max_eigenvalue = eigenvalues[i];
        }
      }
    }
  } else {
    for (uint i = 0; i < numEVexamine; i++) {
      if (eigenvalues[i].imag() > max_eigenvalue.imag()) {
        max_eigenvalue = eigenvalues[i];
      }
    }
  }
  // if (config.branch == BRANCH_TEMPORAL) {
  //   if (config.problem == PB_BOUNDARY_LAYER || config.problem ==
  //   PB_CUSTOM) {
  //     // find the twoeigenvalues with heightes imag part that around
  //     real()=1 double eps = 0.01; bool firstFound = false; bool
  //     secondFound = false; complex maxEV_line, second_maxEV_line; for
  //     (const auto &lambda : eigenvalues) {
  //       if (abs(lambda.real() - 1) < eps && !firstFound) {
  //         maxEV_line = lambda;
  //         firstFound = true;
  //       } else if (abs(lambda.real() - 1) < eps && !secondFound) {
  //         second_maxEV_line = lambda;
  //         secondFound = true;
  //       } else if (abs(lambda.real() - 1) < eps && firstFound &&
  //       secondFound) {
  //         if (lambda.imag() > second_maxEV_line.imag()) {
  //           if (maxEV_line.imag() > lambda.imag()) {
  //             second_maxEV_line = lambda;
  //           } else {
  //             second_maxEV_line = maxEV_line;
  //             maxEV_line = lambda;
  //           }
  //         }
  //       }
  //     }
  //     double min_distance = MAX(maxEV_line.real(),
  //     second_maxEV_line.real()) /
  //                           MIN(maxEV_line.real(),
  //                           second_maxEV_line.real());
  //     std::cout << "maxEV_line = " << maxEV_line.real() << " + "
  //               << maxEV_line.imag() << "i" << std::endl;
  //     std::cout << "second_maxEV_line = " << second_maxEV_line.real() <<
  //     " + "
  //               << second_maxEV_line.imag() << "i" << std::endl;
  //     printf("min_distance = %f\n", min_distance);

  //     double tol = 1 + 1000 * (min_distance - 1);
  //     printf("tol = %f\n", tol);
  //     for (const auto &lambda : eigenvalues) {
  //       // lamb /= config.alpha;
  //       std::cout << "lambda = " << lambda.real() << " + " <<
  //       lambda.imag()
  //                 << "i" << "lambda.real()/max_eigenvalue.real() = " <<
  //                 lambda.real() / max_eigenvalue.real() << std::endl;
  //       if (lambda.real() / max_eigenvalue.real() > tol)
  //         continue; // Skip if real part is greater than 0.75
  //       std::cout << "lambda = " << lambda.real() << " + " <<
  //       lambda.imag()
  //                 << "i" << std::endl;
  //       if (lambda.imag() > max_eigenvalue.imag()) {
  //         max_eigenvalue = lambda;
  //       }
  //     }
  //     std::cout << "max_eigenvalue = " << max_eigenvalue.real() << " + "
  //               << max_eigenvalue.imag() << "i" << std::endl;
  //   }
  // } else {
  //   for (const auto &lambda : eigenvalues) {
  //     if (lambda.imag() > max_eigenvalue.imag()) {
  //       max_eigenvalue = lambda;
  //     }
  //   }
  // }

  // max_eigenvalue = rescaleEV(config.var, max_eigenvalue);
  // std::cout << "max_eigenvalue = " << max_eigenvalue.real() << " + "
  //           << max_eigenvalue.imag() << "i" << std::endl;

  return max_eigenvalue;
}

complex PostProcess::getMostUnstableEigenvalue() const {
  complex max_eigenvalue =
      rescaleEV(config.var, getMostUnstableEigenvalueNotScaled());
  // std::cout << "max_eigenvalue = " << max_eigenvalue.real() << " + "
  //           << max_eigenvalue.imag() << "i" << std::endl;

  return max_eigenvalue;
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
void PostProcess::printSpectrum() const {
  int count = 0;
  int LIMIT = 20; // Limit output to first 20 eigenvalues
  std::string evLabel = config.getEVlabel();
  for (const auto &lambda : eigenvalues) {
    if (count < LIMIT) { // Limit output to first 20 eigenvalues
      auto lambda_scaled = rescaleEV(config.var, lambda);
      std::cout << evLabel << " = " << lambda_scaled.real() << " + "
                << lambda_scaled.imag() << "i" << std::endl;
    }
    count++;
  }
}

void PostProcess::writeToFile(const OSSolver &solver) const {
  // eigenvalues
  std::ofstream file(config.filenameEigenvalues);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << config.filenameEigenvalues
              << std::endl;
    return;
  }
  std::string evLabel = config.getEVlabel();
  file << "# Re(" << evLabel << ")   Im(" << evLabel << ")" << std::endl;
  for (const auto &lambda : eigenvalues) {
    auto lambda_scaled = rescaleEV(config.var, lambda);
    file << lambda_scaled.real() << " " << lambda_scaled.imag() << std::endl;
  }
  file.close();

  // eigenvector
  std::ofstream file2(config.filenameEigenvector);
  if (!file2.is_open()) {
    std::cerr << "Error opening file: " << config.filenameEigenvector
              << std::endl;
    return;
  }
  Eigen::VectorXcd eigenvector =
      solver.computeEigenvector(getMostUnstableEigenvector());

  file2 << "# y   Re(v)   Im(v)" << std::endl;
  for (uint i = 0; i < solver.numQuadPoints; i++) {
    file2 << solver.getYPhysicalRegion(solver, i) << " "
          << eigenvector[i].real() << " " << eigenvector[i].imag() << std::endl;
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
    const auto &lambda = eigenvalues[i];
    auto lambda_scaled = rescaleEV(config.var, lambda);
    const auto &vars = var_complex[i];
    file << lambda_scaled.real() << " " << lambda_scaled.imag() << " "
         << vars.real() << " " << vars.imag() << std::endl;
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
