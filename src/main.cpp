#include "../include/config.hpp"
#include "../include/pp.hpp"
#include "../include/solver.hpp"
#include <chrono>
#include <cstdio>
#include <iostream>

int main() {
  Config config;
  // Load configuration from TOML file
  if (!config.load("config/input.toml")) {
    std::cerr << "Failed to load configuration file." << std::endl;
    return 1;
  }

  std::cout << "Setting up the solver..." << std::endl;
  auto start0 = std::chrono::high_resolution_clock::now();
  OSSolver solver(config);
  auto end0 = std::chrono::high_resolution_clock::now();
  printf("Setup time: %.2f seconds\n",
         std::chrono::duration<double>(end0 - start0).count());

  std::vector<complex> eigenvalues;
  EigenSolution eig;

  std::string varLabel = config.getVarlabel(config.branch);

  if (config.multipleRun) {
    std::vector<complex> vars;
    complex evmax;
    double dvar_r =
        (config.vars_r.num == 1)
            ? 0
            : (config.vars_r.max - config.vars_r.min) / (config.vars_r.num - 1);
    double dvar_i =
        (config.vars_i.num == 1)
            ? 0
            : (config.vars_i.max - config.vars_i.min) / (config.vars_i.num - 1);
    for (int i = 0; i < config.vars_r.num; i++) {
      for (int j = 0; j < config.vars_i.num; j++) {
        double var_r = config.vars_r.min + i * dvar_r;
        double var_i = config.vars_i.min + j * dvar_i;
        complex var(var_r, var_i);
        std::cout << "Running simulation for " << varLabel << " = " << var
                  << ", β = " << config.beta << ", Re = " << config.re;
        config.setVar(var);
        solver.buildMatrices(config.branch);
        eig = solver.solve();
        PostProcess pp(eig, config);
        evmax = pp.getTargetEVal();

        std::cout << " ==> Most unstable EV: " << evmax << std::endl;

        eigenvalues.push_back(evmax);
        vars.push_back(var);
      }
    }
    // Print the results
    PostProcess::writeToFile(vars, eigenvalues, config);
  } else {
    std::cout << "Running simulation for " << varLabel << " = " << config.var
              << ", β = " << config.beta << ", Re = " << config.re << std::endl;

    {
      std::cout << "Building matrices..." << std::endl;
      auto start = std::chrono::high_resolution_clock::now();
      solver.buildMatrices(config.branch);
      auto end = std::chrono::high_resolution_clock::now();
      printf("Build time: %.2f seconds\n",
             std::chrono::duration<double>(end - start).count());
    }
    {
      std::cout << "Running simulation..." << std::endl;
      auto start = std::chrono::high_resolution_clock::now();
      eig = solver.solve();
      auto end = std::chrono::high_resolution_clock::now();
      printf("Solve time: %.2f seconds\n",
             std::chrono::duration<double>(end - start).count());
    }
    PostProcess pp(eig, config);

    // Print the results
    pp.printSpectrum();

    // Print the most unstable eigenvalue
    complex mostUnstableEV = pp.getTargetEVal();

    // more decimals
    std::cout.precision(15);
    std::cout << "Most unstable eigenvalue: " << mostUnstableEV << std::endl;
    // std::cout << "Most unstable eigenvalue: " << mostUnstableEV.real() << " +
    // "
    //           << mostUnstableEV.imag() << "i" << std::endl;

    Vector evec = solver.computeEVec(pp.getTargetEVec(), config.branch);

    pp.writeToFile(solver.gaussPoints, evec, solver.Uprof);
  }

  return 0;
}
