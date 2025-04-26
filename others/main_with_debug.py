from config import Config, ProblemType, Branch
import os
from methods import runUFromFile, runPoiseuille
import numpy as np
import matplotlib.pyplot as plt
from OS import Orr_Sommerfeld_Espacial, Orr_Sommerfeld_Temporal


def main():
    # read data
    CONFIG_FILE = "data/input.toml"
    DIR_SCRIPT = os.path.dirname(os.path.realpath(__file__))
    CONFIG_FILE = os.path.join(DIR_SCRIPT, "../", CONFIG_FILE)

    config = Config.from_toml(CONFIG_FILE)

    if config.run_multiple:
        max_eigs = []
        for i, var_i in enumerate(config.vars_range_i):
            aux = []
            for j, var_r in enumerate(config.vars_range_r):
                print(
                    f"real {j + 1}/{len(config.vars_range_r)} imag {i + 1}/{len(config.vars_range_i)}"
                )
                maxeig = 0
                if config.branch == Branch.Temporal:
                    config.alpha = complex(var_r, var_i)
                else:
                    config.omega = complex(var_r, var_i)

                if config.problem == ProblemType.Poiseuille:
                    maxeig = runPoiseuille(config)
                else:
                    maxeig = runUFromFile(config)
                aux.append(maxeig)
            max_eigs.append(aux)

        # plot curves of points (real and imag part of the points) of the same color for same alpha_i/different alpha_r
        max_eigs = np.array(max_eigs)

        for i, vars_i in enumerate(config.vars_range_i):
            eigs = max_eigs[i]  # Eigenvalues corresponding to the same imaginary part
            real_parts = np.real(eigs)
            imag_parts = np.imag(eigs)
            plt.scatter(
                real_parts,
                imag_parts,
                label=f"Im(α) = {np.round(vars_i[i], 2)}",
                alpha=0.7,
            )

        print("alphas_r", config.vars_range_r)
        print("alphas_i", config.vars_range_i)
        plt.xlabel(config.plot_label + "_r")
        plt.ylabel(config.plot_label + "_i")
        plt.legend()
        plt.grid()
        plt.show()
    else:
        if config.problem == ProblemType.Poiseuille:
            runPoiseuille(config)
            if config.debug:
                print("Reference code:")
                if config.branch == Branch.Temporal:
                    Orr_Sommerfeld_Temporal(
                        config.n - 1, config.re, config.alpha, config.beta, 0, 0
                    )
                else:
                    Orr_Sommerfeld_Espacial(
                        config.n - 1, config.re, config.omega, config.beta, 0, 0
                    )
        else:
            runUFromFile(config)
            if config.debug:
                print("Reference code:")
                if config.branch == Branch.Temporal:
                    Orr_Sommerfeld_Temporal(
                        config.n - 1, config.re, config.alpha, config.beta, 0, 2
                    )
                else:
                    Orr_Sommerfeld_Espacial(
                        config.n - 1, config.re, config.omega, config.beta, 0, 2
                    )


if __name__ == "__main__":
    main()
