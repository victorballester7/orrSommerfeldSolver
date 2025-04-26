#include "../include/orr_sommerfeld.hpp"
#include "../include/pp.hpp"
#include "../include/read_conf.hpp"
#include <chrono>
#include <cstdio>
#include <iostream>

std::vector<complex> runSim(OSSolver &solver, bool print = false) {
  // Create and set up the solver
  std::chrono::high_resolution_clock::time_point start, end;
  if (print) {
    // Solve the generalized eigenvalue problem
    std::cout << "Solving the eigenvalue problem..." << std::endl;
    start = std::chrono::high_resolution_clock::now();
  }
  std::vector<complex> eigenvalues = solver.solve();
  if (print) {
    end = std::chrono::high_resolution_clock::now();
    printf("Solve time: %.2f seconds\n",
           std::chrono::duration<double>(end - start).count());
    // Placeholder for the actual simulation logic
  }
  return eigenvalues;
}

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

  if (config.run_multiple) {
    std::vector<complex> ev_aux, vars;
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
        complex var_print = (config.problem == PB_BOUNDARY_LAYER)
                                ? (var * DELTASTAR_BLASIUS)
                                : var;
        std::cout << "Running simulation for var = " << var_print << std::endl;
        config.setVar(var);
        solver.setVar(var, config.branch);
        solver.buildMatrices(config.branch);
        ev_aux = runSim(solver);
        PostProcess pp(config, ev_aux);
        evmax = pp.getMostUnstableEigenvalue();
        // std::cout << "Most unstable eigenvalue: " << evmax.real()
        //       << " + " << evmax.imag() << "i" << std::endl;

        eigenvalues.push_back(evmax);
        vars.push_back(var_print);
      }
    }
    // Print the results
    PostProcess pp(config, eigenvalues, false);
    // pp.plotSpectrum();
    pp.writeToFile(config.filenameEigenvalues, vars);

  } else {
    std::cout << "Running simulation for "
              << ((config.branch == BRANCH_TEMPORAL) ? "α = " : "ω = ")
              << config.var << std::endl;
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
      eigenvalues = runSim(solver);
      auto end = std::chrono::high_resolution_clock::now();
      printf("Solve time: %.2f seconds\n",
             std::chrono::duration<double>(end - start).count());
    }
    PostProcess pp(config, eigenvalues);

    // Print the results
    // pp.printSpectrum(eigenvalues);
    pp.printSpectrum();

    // Print the most unstable eigenvalue
    complex mostUnstableEigenvalue = pp.getMostUnstableEigenvalue();
    std::cout << "Most unstable eigenvalue: " << mostUnstableEigenvalue.real()
              << " + " << mostUnstableEigenvalue.imag() << "i" << std::endl;

    pp.writeToFile(config.filenameEigenvalues);
    // if (config.doPlot) {
    //   pp.plotSpectrum();
    // }
  }

  return 0;
}
