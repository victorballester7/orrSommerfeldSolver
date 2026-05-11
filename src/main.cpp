#include "../include/config.hpp"
#include "../include/pp.hpp"
#include "../include/solver.hpp"
#include <chrono>
#include <cstdio>
#include <iomanip>
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
    std::string varLabel = config.getVarlabel();
    for (int i = 0; i < config.vars_r.num; i++) {
      for (int j = 0; j < config.vars_i.num; j++) {
        double var_r = config.vars_r.min + i * dvar_r;
        double var_i = config.vars_i.min + j * dvar_i;
        complex var(var_r, var_i);
        std::cout << std::fixed << std::setprecision(6)
                  << "Running simulation for " << varLabel << " = " << var
                  << ", β = " << config.beta << ", Re = " << config.re
                  << std::endl;
        config.setVar(var);
        std::cout << "config.var: " << config.var << std::endl;
        std::cout << "config.k2: " << config.k2 << std::endl;
        std::cout << "solver.var: " << solver.var << std::endl;
        std::cout << "solver.k2: " << solver.k2 << std::endl;
        solver.setVar(var, config.branch);
        std::cout << "solver.var: " << solver.var << std::endl;
        std::cout << "solver.k2: " << solver.k2 << std::endl;
        solver.buildMatrices(config.branch);
        eig = solver.solve();
        PostProcess pp(eig);
        // pp.printSpectrum();
        // evmax = pp.getMostUnstableEigenvalueNotScaled();
        evmax = pp.getMostUnstableEigenvalue();
        // std::cout << "Most unstable eigenvalue: " << evmax.real()
        //       << " + " << evmax.imag() << "i" << std::endl;

        eigenvalues.push_back(evmax);
        vars.push_back(var);
      }
    }
    // Print the results
    PostProcess pp(eigenvalues);
    // pp.plotSpectrum();
    pp.writeToFile(vars);

  } else {
    std::cout << std::fixed << std::setprecision(6) << "Running simulation for "
              << config.getVarlabel() << " = " << config.var
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
    PostProcess pp(eig);

    std::string evLabel = config.getEVlabel();

    // Print the results
    pp.printSpectrum(evLabel);

    // Print the most unstable eigenvalue
    complex mostUnstableEigenvalue = pp.getMostUnstableEigenvalue();
    // more decimals
    std::cout.precision(15);
    std::cout << "Most unstable eigenvalue: " << mostUnstableEigenvalue.real()
              << " + " << mostUnstableEigenvalue.imag() << "i" << std::endl;

    pp.writeToFile(evLabel, solver);
  }

  return 0;
}
