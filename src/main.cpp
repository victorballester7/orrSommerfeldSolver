#include "../include/config.hpp"
#include "../include/pp.hpp"
#include "../include/solver.hpp"
#include <chrono>
#include <cstdio>
#include <iostream>

int main() {
  // Parameters - can be modified as needed
  // int p = 50;       // Polynomial degree (e.g., 100, 200, 500)
  // double Re = 2000; // Reynolds number
  // double alpha = 1;  // Wavenumber

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
  Eigen::ComplexEigenSolver<Matrix> eig;

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
        complex var_print = (config.problem == PB_BOUNDARY_LAYER)
                                ? (var * DELTASTAR_BLASIUS)
                                : var;
        complex beta_print = (config.problem == PB_BOUNDARY_LAYER)
                                  ? (config.beta * DELTASTAR_BLASIUS)
                                  : config.beta;
        std::cout << "Running simulation for " << varLabel << " = " << var_print
                  << ", β = " << beta_print << std::endl;
        config.setVar(var);
        solver.setVar(var, config.branch);
        solver.buildMatrices(config.branch);
        eig = solver.solve();
        PostProcess pp(config, eig);
        // pp.printSpectrum();
        evmax = pp.getMostUnstableEigenvalueNotScaled();
        // std::cout << "Most unstable eigenvalue: " << evmax.real()
        //       << " + " << evmax.imag() << "i" << std::endl;

        eigenvalues.push_back(evmax);
        vars.push_back(var_print);
      }
    }
    // Print the results
    PostProcess pp(config, eigenvalues);
    // pp.plotSpectrum();
    pp.writeToFile(vars);

  } else {
    complex var_print = (config.problem == PB_BOUNDARY_LAYER)
                            ? (config.var * DELTASTAR_BLASIUS)
                            : config.var;
    complex beta_print = (config.problem == PB_BOUNDARY_LAYER)
                            ? (config.beta * DELTASTAR_BLASIUS)
                            : config.beta;

    std::cout << "Running simulation for " << config.getVarlabel() << " = "
              << var_print << ", β = " << beta_print << std::endl;
    std::cout << "Building matrices..." << std::endl;

    {
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
    PostProcess pp(config, eig);

    // Print the results
    // pp.printSpectrum(eigenvalues);
    pp.printSpectrum();

    // Print the most unstable eigenvalue
    complex mostUnstableEigenvalue = pp.getMostUnstableEigenvalue();
    // more decimals
    std::cout.precision(15);
    std::cout << "Most unstable eigenvalue: " << mostUnstableEigenvalue.real()
              << " + " << mostUnstableEigenvalue.imag() << "i" << std::endl;

    pp.writeToFile(solver);
    // if (config.doPlot) {
    //   pp.plotSpectrum();
    // }
  }

  return 0;
}
