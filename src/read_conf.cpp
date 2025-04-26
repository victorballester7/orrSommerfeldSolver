#include "../include/read_conf.hpp"
#include "../libs/toml.hpp" // Make sure toml.hpp is in your include path
#include <iostream>
#include <string>

// Load config from TOML file
bool Config::load(const std::string &filename) {
  try {
    auto tbl = toml::parse_file(filename);

    // General
    const auto &general = tbl["general"];
    p = (uint)general["n"].value_or(0);
    re = general["re"].value_or(0.0);
    alpha = parseComplex(general["alpha"]);
    beta = parseComplex(general["beta"]);
    omega = parseComplex(general["omega"]);

    // Flags
    const auto &flags = tbl["flags"];
    branch = flags["branch"].value_or("");
    if (!isValid(branch, branches)) {
      std::cerr << "Invalid branch: " << branch
        << ". Available options are: ";
      for (const auto &b : branches) {
        std::cerr << b << " ";
      }
      std::cerr << std::endl;
      return false;
    }
    problem = flags["problem"].value_or("");
    if (!isValid(problem, problems)) {
      std::cerr << "Invalid problem: " << problem
        << ". Available options are: ";
      for (const auto &_p : problems) {
        std::cerr << _p << " ";
      }
      std::cerr << std::endl;
      return false;
    }
    filenameEigenvalues = flags["fileWriteEigenvalues"].value_or("");
    doPlot = flags["doPlot"].value_or(false);
    use_c = flags["use_c"].value_or(false);
    run_multiple = flags["run_multiple"].value_or(false);


    // Custom problem flags
    const auto &customProblemFlags = tbl["customProblemFlags"];
    filenameUprofile = customProblemFlags["filenameUprofile"].value_or("");
    colX = static_cast<uint>(customProblemFlags["colX"].value_or(0));
    colY = static_cast<uint>(customProblemFlags["colY"].value_or(0));
    numSkipHeaderLines = static_cast<uint>(customProblemFlags["numSkipHeaderLines"].value_or(0));


    // Vars
    const auto &runMultipleFlags = tbl["runMultipleFlags"];
    vars_r = parseRange(runMultipleFlags["vars_r"]);
    vars_i = parseRange(runMultipleFlags["vars_i"]);

    // Plot
    const auto &plot = tbl["plot"];
    plot_lims = parsePlotLims(plot["plot_lims"]);

    // rescale quantities for Blasius flow only
    deltaStarRescaling();

    k2 = getK2();
    

    return true;
  } catch (const toml::parse_error &err) {
    std::cerr << "TOML Parse error: " << err << std::endl;
    return false;
  }
}
