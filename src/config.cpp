#include "../include/config.hpp"
#include "../libs/toml.hpp" // Make sure toml.hpp is in your include path
#include <iostream>
#include <string>

// Load config from TOML file
bool Config::load(const std::string &filename) {
  try {
    auto tbl = toml::parse_file(filename);

    // Params
    const auto &params = tbl["params"];
    p = (uint)params["n"].value_or(0);
    re = params["re"].value_or(0.0);
    beta = parseComplex(params["beta"]);

    // single run
    const auto &singleRunParams = tbl["singleRunParams"];
    var = parseComplex(singleRunParams["var"]);

    // multiple run
    const auto &multipleRunParams = tbl["multipleRunParams"];
    vars_r = parseRange(multipleRunParams["vars_r"]);
    vars_i = parseRange(multipleRunParams["vars_i"]);
    
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
    multipleRun = flags["multipleRun"].value_or(false);

    // Custom problem flags
    const auto &customProblemFlags = tbl["customProblemFlags"];
    filenameUprofile = customProblemFlags["filenameUprofile"].value_or("");
    colX = static_cast<uint>(customProblemFlags["colX"].value_or(0));
    colY = static_cast<uint>(customProblemFlags["colY"].value_or(0));
    numSkipHeaderLines =
        static_cast<uint>(customProblemFlags["numSkipHeaderLines"].value_or(0));
    // Plot
    const auto &plot = tbl["plot"];
    plotLims = parsePlotLims(plot["plotLims"]);

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
