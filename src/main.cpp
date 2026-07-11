#include "../include/config.hpp"
#include "../include/pp.hpp"
#include "../include/solver.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <exception>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

int main() {
  Config config;
  // Load configuration from TOML file
  if (!config.load("config/input.toml")) {
    std::cerr << "Failed to load configuration file." << std::endl;
    return 1;
  }

  std::vector<complex> eigenvalues;
  EigenSolution eig;

  std::string varLabel = config.getVarlabel(config.branch);

  if (config.multipleRun) {
    double dvar_r =
        (config.vars_r.num == 1)
            ? 0
            : (config.vars_r.max - config.vars_r.min) / (config.vars_r.num - 1);
    double dvar_i =
        (config.vars_i.num == 1)
            ? 0
            : (config.vars_i.max - config.vars_i.min) / (config.vars_i.num - 1);

    const int numVarsR = config.vars_r.num;
    const int numVarsI = config.vars_i.num;
    const int totalRuns = numVarsR * numVarsI;
    const unsigned int hardwareThreads = std::thread::hardware_concurrency();
    const int availableThreads =
        static_cast<int>(hardwareThreads ? hardwareThreads : 1);
    const int numThreads = std::max(1, std::min(totalRuns, availableThreads));

    std::vector<complex> vars(totalRuns);
    eigenvalues.resize(totalRuns);

    std::atomic<int> nextRun{0};
    std::atomic<bool> stopRequested{false};
    std::exception_ptr firstException = nullptr;
    std::mutex exceptionMutex;
    std::mutex coutMutex;
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(numThreads));

    auto runPoint = [&](int runIdx, OSSolver &workerSolver,
                        Config &workerConfig) {
      const int i = runIdx / numVarsI;
      const int j = runIdx % numVarsI;
      double var_r = config.vars_r.min + i * dvar_r;
      double var_i = config.vars_i.min + j * dvar_i;
      complex var(var_r, var_i);

      workerConfig.setVar(var);
      workerSolver.buildMatrices(workerConfig.branch);
      EigenSolution workerEig = workerSolver.solve();
      PostProcess pp(workerEig, workerConfig);
      complex evmax = pp.getTargetEVal();

      {
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << "Running simulation for " << varLabel << " = " << var
                  << ", β = " << workerConfig.beta << ", Re = "
                  << workerConfig.re << " ==> Most unstable EV: " << evmax
                  << std::endl;
      }

      vars[runIdx] = var;
      eigenvalues[runIdx] = evmax;
    };

    auto worker = [&]() {
      Config workerConfig(config);
      OSSolver workerSolver(workerConfig);

      while (true) {
        if (stopRequested) {
          return;
        }

        const int runIdx = nextRun.fetch_add(1);
        if (runIdx >= totalRuns) {
          return;
        }

        try {
          runPoint(runIdx, workerSolver, workerConfig);
        } catch (...) {
          std::lock_guard<std::mutex> lock(exceptionMutex);
          if (!firstException) {
            firstException = std::current_exception();
            stopRequested = true;
          }
          return;
        }
      }
    };

    std::cout << "Running " << totalRuns << " simulations using " << numThreads
              << " worker thread(s)..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    for (int threadIdx = 0; threadIdx < numThreads; ++threadIdx) {
      workers.emplace_back(worker);
    }

    for (auto &thread : workers) {
      thread.join();
    }

    if (firstException) {
      std::rethrow_exception(firstException);
    }

    auto end = std::chrono::high_resolution_clock::now();
    printf("Multiple run time: %.2f seconds\n",
           std::chrono::duration<double>(end - start).count());

    // Print the results
    PostProcess::writeToFile(vars, eigenvalues, config);
  } else {
    std::cout << "Setting up the solver..." << std::endl;
    auto start0 = std::chrono::high_resolution_clock::now();
    OSSolver solver(config);
    auto end0 = std::chrono::high_resolution_clock::now();
    printf("Setup time: %.2f seconds\n",
           std::chrono::duration<double>(end0 - start0).count());

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
