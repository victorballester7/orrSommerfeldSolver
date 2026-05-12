#include "../include/config.hpp"
#include "../libs/toml.hpp" // Make sure toml.hpp is in your include path
#include <iostream>
#include <string>

std::ostream &operator<<(std::ostream &os, const complex &z) {
  const double re = z.real();
  const double im = z.imag();
  const double tol = 1e-10;

  if (std::abs(re) > tol) {
    os << re;
  }

  if (im > tol) {
    os << " + " << im << "i";
  } else if (im < -tol) {
    os << " - " << std::abs(im) << "i";
  } 

  if (std::abs(re) <= tol && std::abs(im) <= tol) {
    os << "0";
  }

  return os;
}

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
    targetEV = parseComplex(singleRunParams["targetEV"]);
    useTargetEV = singleRunParams["useTargetEV"].value_or(false);

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

    doPlot = flags["doPlot"].value_or(false);
    use_c = flags["use_c"].value_or(false);
    multipleRun = flags["multipleRun"].value_or(false);

    // Custom problem flags
    const auto &customProblemFlags = tbl["customProblemFlags"];
    filenameUprofile = customProblemFlags["filenameUprofile"].value_or("");
    blasiusLikeDomain = customProblemFlags["blasiusLikeDomain"].value_or(false);
    colX = static_cast<uint>(customProblemFlags["colX"].value_or(0));
    colY = static_cast<uint>(customProblemFlags["colY"].value_or(0));
    numSkipHeaderLines =
        static_cast<uint>(customProblemFlags["numSkipHeaderLines"].value_or(0));

    // Plot
    const auto &plot = tbl["plot"];
    plotLims = parsePlotLims(plot["plotLims"]);

    k2 = getK2();

    return true;
  } catch (const toml::parse_error &err) {
    std::cerr << "TOML Parse error: " << err << std::endl;
    return false;
  }
}
