#include "../include/config.hpp"
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
    var = parseComplex(general["var"]);
    beta = parseComplex(general["beta"]);

    // Flags
    const auto &flags = tbl["flags"];
    std::string field = "branch";

    branch = flags[field].value_or("");
    isValidMessage(field, branch, branches);

    field = "problem";
    problem = flags[field].value_or("");
    isValidMessage(field, problem, problems);

    filenameEigenvalues = flags["fileWriteEigenvalues"].value_or("");
    filenameEigenvector = flags["fileWriteEigenvector"].value_or("");
    doPlot = flags["doPlot"].value_or(false);
    use_c = flags["use_c"].value_or(false);
    run_multiple = flags["run_multiple"].value_or(false);

    // Custom problem flags
    const auto &customProblemFlags = tbl["customProblemFlags"];
    filenameUprofile = customProblemFlags["filenameUprofile"].value_or("");
    colX = static_cast<uint>(customProblemFlags["colX"].value_or(0));
    colY = static_cast<uint>(customProblemFlags["colY"].value_or(0));
    numSkipHeaderLines =
        static_cast<uint>(customProblemFlags["numSkipHeaderLines"].value_or(0));

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

std::string Config::getEVlabel() const {
  if (branch == BRANCH_TEMPORAL) {
    if (use_c)
      return "c";
    else
      return "ω";
  } else {
    return "α";
  }
}

std::string Config::getVarlabel() const {
  if (branch == BRANCH_TEMPORAL) {
    return "α";
  } else {
    return "ω";
  }
}
