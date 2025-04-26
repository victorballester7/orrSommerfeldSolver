#include "../include/pp.hpp"
#include <matplot/matplot.h>

void PostProcess::rescaleEVector(const std::complex<double> &alpha) {
  // Rescale eigenvalues by the factor alpha
  for (auto &lambda : eigenvalues) {
    lambda /= alpha;
  }
}

complex PostProcess::getMostUnstableEigenvalue() {
  // Find the eigenvalue with the largest imaginary part
  complex max_eigenvalue = eigenvalues[0];
  for (const auto &lambda : eigenvalues) {
    if (config.branch == BRANCH_TEMPORAL) {
      if (config.problem == PB_BOUNDARY_LAYER || config.problem == PB_CUSTOM) {
        // we work with the rescaled eigenvalues
        // lamb /= config.alpha;
        if (lambda.real() / max_eigenvalue.real() > 2.0)
          continue; // Skip if real part is greater than 0.75
      }
    }
    if (lambda.imag() > max_eigenvalue.imag()) {
      max_eigenvalue = lambda;
    }
  }

  return max_eigenvalue;
}

// Print the spectrum
void PostProcess::printSpectrum() {
  int count = 0;
  int LIMIT = 20; // Limit output to first 20 eigenvalues
  for (const auto &lambda : eigenvalues) {
    if (count < LIMIT) { // Limit output to first 20 eigenvalues
      if (config.branch == BRANCH_TEMPORAL) {
        if (config.use_c) {
          std::cout << "c = ";
        } else {
          std::cout << "ω = ";
        }
      } else {
        std::cout << "α = ";
      }
      std::cout << lambda.real() << " + " << lambda.imag() << "i" << std::endl;
    }
    count++;
  }
}

void PostProcess::writeToFile(const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (const auto &lambda : eigenvalues) {
    file << lambda.real() << " " << lambda.imag() << std::endl;
  }

  file.close();
}

void PostProcess::writeToFile(const std::string &filename,
                              const std::vector<complex> &var_complex) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (uint i = 0; i < eigenvalues.size(); i++) {
    const auto &lambda = eigenvalues[i];
    const auto &vars = var_complex[i];
    file << lambda.real() << " " << lambda.imag() << " " << vars.real() << " "
         << vars.imag() << std::endl;
  }

  file.close();
}

void PostProcess::plotSpectrum() {
  using namespace matplot;
  // auto theta = linspace(0, 1, 500);
  // auto x = transform(
  //     theta, [&](double theta) { return exp(theta) * sin(100 * theta); });
  // auto y = transform(
  //     theta, [&](double theta) { return exp(theta) * cos(100 * theta); });

  // auto s = scatter(x, y);
  // s->marker_color("b");
  // s->marker_face_color({0, .5, .5});
  // At the beginning of your function
  // show();
  // Extract real and imaginary parts
  std::vector<double> re, im;
  for (const auto &eig : eigenvalues) {
    re.push_back(eig.real()); // Real part
    im.push_back(eig.imag()); // Imaginary part
  }

  // Create a scatter plot
  // for (int i = 0; i < re.size(); i += chunk_size_points) {
  //   auto re_aux = std::vector<double>(re.begin() + i,
  //                         re.begin() + std::min(i + chunk_size_points,
  //                         (int)re.size()));
  //   auto im_aux = std::vector<double>(im.begin() + i,
  //                         im.begin() + std::min(i + chunk_size_points,
  //                         (int)im.size()));
  //   auto s = scatter(re_aux, im_aux);
  // }

  auto s = scatter(re, im); // Plot with red color and size 10
  s->marker_color("b");
  s->marker_face_color({0, .5, .5});

  // add grid
  grid(true);

  // Set axis limits
  xlim({config.plot_lims.xmin, config.plot_lims.xmax});
  ylim({config.plot_lims.ymin, config.plot_lims.ymax});

  show();
}
