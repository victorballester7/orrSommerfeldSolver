#ifndef READ_CONF_HPP
#define READ_CONF_HPP

#include "../libs/toml.hpp" // Make sure toml.hpp is in your include path
#include <complex>
#include <iostream>
#include <set>
#include <string>

using complex = std::complex<double>;

// define available values for branch: temporal or spatial

constexpr const char *BRANCH_TEMPORAL = "temporal";
constexpr const char *BRANCH_SPATIAL = "spatial";
constexpr const char *PB_POISEUILLE = "Poiseuille";
constexpr const char *PB_BOUNDARY_LAYER = "BoundaryLayer";
constexpr const char *PB_COUETTE = "Couette";
constexpr const char *PB_CUSTOM = "Custom";

#define DELTASTAR_BLASIUS 1.7207876573

struct Range {
  double min;
  double max;
  int num;
};

struct PlotLims {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

class Config {
public:
  // General
  uint p;
  double re;
  complex var;
  complex beta;
  complex k2;

  // Flags
  std::string branch;
  const std::set<std::string> branches = {BRANCH_TEMPORAL, BRANCH_SPATIAL};
  std::string problem;
  const std::set<std::string> problems = {PB_POISEUILLE, PB_BOUNDARY_LAYER,
                                          PB_COUETTE, PB_CUSTOM};
  std::string filenameEigenvalues;
  std::string filenameEigenvector;
  bool doPlot;
  bool use_c;
  bool run_multiple;

  // Custom problem flags
  std::string filenameUprofile;
  uint colX, colY;
  uint numSkipHeaderLines;

  // Run multiple flags
  Range vars_r;
  Range vars_i;

  // Plot
  PlotLims plot_lims;

  // Load config from TOML file
  bool load(const std::string &filename);

  // Set the variable for the simulation
  void setVar(complex var_) {
    var = var_;
    if (branch == BRANCH_TEMPORAL) {
      k2 = var * var + beta * beta;
    }
  }

  std::string getVarlabel() const;
  std::string getEVlabel() const;

private:
  void deltaStarRescaling() {
    if (problem == PB_BOUNDARY_LAYER) {
      // rescale quantities for Blasius flow because delta* is not 1, it is
      // DELTASTAR_BLASIUS
      var = var / DELTASTAR_BLASIUS;
      beta = beta / DELTASTAR_BLASIUS;

      vars_r.min = vars_r.min / DELTASTAR_BLASIUS;
      vars_r.max = vars_r.max / DELTASTAR_BLASIUS;
      vars_i.min = vars_i.min / DELTASTAR_BLASIUS;
      vars_i.max = vars_i.max / DELTASTAR_BLASIUS;

      re = re / DELTASTAR_BLASIUS;
    }
  }

  complex parseComplex(const toml::node_view<toml::node> &node) const {
    complex c;
    c.real(node["r"].value_or(0.0));
    c.imag(node["i"].value_or(0.0));
    return c;
  }

  Range parseRange(const toml::node_view<toml::node> &node) const {
    Range r;
    r.min = node["min"].value_or(0.0);
    r.max = node["max"].value_or(0.0);
    r.num = node["num"].value_or(1);
    if (r.num == 0) {
      r.num = 1;
    }
    return r;
  }

  PlotLims parsePlotLims(const toml::node_view<toml::node> &node) const {
    PlotLims plt;
    plt.xmin = node["xmin"].value_or(0.0);
    plt.xmax = node["xmax"].value_or(0.0);
    plt.ymin = node["ymin"].value_or(0.0);
    plt.ymax = node["ymax"].value_or(0.0);
    return plt;
  }

  complex getK2() const {
    // Calculate k2 based on alpha and beta
    return std::pow(var, 2) + std::pow(beta, 2);
  }

  bool isValid(const std::string &value,
               const std::set<std::string> &validValues) const {
    return validValues.find(value) != validValues.end();
  }

  bool isValidMessage(const std::string field, const std::string &value,
                      const std::set<std::string> &validValues) const {
    if (!isValid(value, validValues)) {
      std::cerr << "Invalid " << field << ": " << value
                << ". Available options are: ";
      for (const auto &v : validValues) {
        std::cerr << v << " ";
      }
      std::cerr << std::endl;
      return false;
    }
    return true;
  }
};

#endif // READ_CONF_HPP
