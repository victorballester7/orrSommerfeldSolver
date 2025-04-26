#ifndef PP_HPP
#define PP_HPP

#include "read_conf.hpp"
#include <complex>
#include <vector>

class PostProcess {
  Config &config;
  std::vector<complex> &eigenvalues;

  void rescaleEVector(const std::complex<double> &alpha);

public:
  PostProcess(Config &_config, std::vector<complex> &_eigenvalues,
              bool scaling = true)
      : config(_config), eigenvalues(_eigenvalues) {
    if (scaling && config.branch == BRANCH_TEMPORAL && config.use_c && abs(config.alpha) > 1e-10) {
      rescaleEVector(config.alpha);
    }
  };

  complex getMostUnstableEigenvalue();

  // Print the spectrum
  void printSpectrum();

  // Write the spectrum to a file
  void writeToFile(const std::string &filename);
  void writeToFile(const std::string &filename,
                              const std::vector<complex> &vars);
  // Plot the spectrum
  void plotSpectrum();
};
#endif // PP_HPP
