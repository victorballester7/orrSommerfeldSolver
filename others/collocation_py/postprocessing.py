import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from config import PlotLimits, ProblemType, Branch, Config
from typing import Tuple

# def printEVSorted(EVvariable, maxEVs, isBL, isTemporal):
def printEVSorted(EVvariable: npt.NDArray[np.complex128], maxEVs: int, problem: ProblemType, branch: Branch) -> None:
    EVaux = EVvariable
    if branch == Branch.Temporal:
        idx = np.argsort(-EVvariable.imag)
    else:
        # filter for extra large eigenvalues
        EVaux = EVaux[abs(EVaux) < 5]
        idx = np.argsort(-EVaux.real)

    EVaux = EVaux[idx]
    if problem == ProblemType.BoundaryLayer or problem == ProblemType.Gap:
        threshold_realpart = 0.95
        EVaux = EVaux[abs(EVaux.real) < threshold_realpart]

    maxEV, _ = getMostUnstableEV(EVaux, np.zeros((EVaux.size, EVaux.size)).astype(np.complex128), problem, branch)

    print(f"Most unstable EV: {maxEV.real}    {maxEV.imag}i")

    for i in range(min(maxEVs, len(EVaux))):
        print(f"EV {i}: {EVaux[i].real}    {EVaux[i].imag}i")
    print("...")


def updateLimits(problem: ProblemType,branch: Branch, plotLimits: PlotLimits) -> PlotLimits:
    if branch == Branch.Temporal:
        if problem == ProblemType.BoundaryLayer or problem == ProblemType.Gap:
            plotLimits.xmin = 0.1
            plotLimits.xmax = 1.1
            plotLimits.ymin = -1
            plotLimits.ymax = 0.1
        else:
            plotLimits.xmin = 0.1
            plotLimits.xmax = 1
            plotLimits.ymin = -1
            plotLimits.ymax = 0.1
    else:
        if problem == ProblemType.BoundaryLayer or problem == ProblemType.Gap:
            plotLimits.xmin = -0.2
            plotLimits.xmax = 1
            plotLimits.ymin = -0.1
            plotLimits.ymax = 0.8
        else:
            plotLimits.xmin = -0.2
            plotLimits.xmax = 1
            plotLimits.ymin = -0.1
            plotLimits.ymax = 0.1

    return plotLimits



def getMostUnstableEV(EVvariable: npt.NDArray[np.complex128], vv: npt.NDArray[np.complex128], problem: ProblemType, branch: Branch) -> Tuple[complex, npt.NDArray[np.complex128]]:
    EVaux = EVvariable
    if branch == Branch.Temporal:
        EVaux = EVaux[np.abs(EVvariable) < 1.5]
        # if for each EV the minimum distance to another EV is larger than 0.1, ignore it
        for i in range(len(EVaux)):
            if np.min(np.abs(EVaux - EVaux[i])) > 0.2:
                EVaux[i] = -99j
    i = np.argmax(EVaux.imag)
    return EVaux[i], vv[:, i]


def getPerturbation(v: npt.NDArray[np.complex128], eta: npt.NDArray[np.complex128], alpha: complex, beta: complex, D1: npt.NDArray[np.complex128]) -> Tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    k2 = alpha**2 + beta**2

    u = 1.0j / k2 * (alpha * D1 @ v - beta * eta)
    w = 1.0j / k2 * (beta * D1 @ v + alpha * eta)
    return u, v, w


def plotSpectrum(EVvariable: npt.NDArray[np.complex128], EVlabel: str, plotLimits: PlotLimits) -> None:
    # Define and train Kernel Ridge Regression model
    plt.plot(EVvariable.real, EVvariable.imag, "ob", markersize=5)
    
    plt.xlabel(f"{EVlabel}_r")
    plt.ylabel(f"{EVlabel}_i")
    # xmin = np.max(np.array([np.min(EVvariable.real), -10]))
    # xmax = np.min(np.array([np.max(EVvariable.real), 10]))
    # ymin = np.max(np.array([np.min(EVvariable.imag), -0.4]))
    # ymax = np.min(np.array([np.max(EVvariable.imag), 1]))

    plt.xlim([plotLimits.xmin, plotLimits.xmax])
    plt.ylim([plotLimits.ymin, plotLimits.ymax])
    # plt.xlim([-0.2, 1])
    # plt.ylim([-1, 0.1])
    plt.grid()
    plt.show()


def plotEVector(y_cheb: npt.NDArray[np.float64], u: npt.NDArray[np.complex128], v: npt.NDArray[np.complex128]) -> None:
    plt.plot(u.real, y_cheb, label="u_r")
    plt.plot(u.imag, y_cheb, label="u_i")

    plt.plot(v.real, y_cheb, label="v_r")
    plt.plot(v.imag, y_cheb, label="v_i")

    plt.xlabel("v")
    plt.ylabel("y")
    plt.legend()
    plt.show()


def postprocess(
        eigs: npt.NDArray[np.complex128], config: Config) -> None: 
    maxEVs = 10
    printEVSorted(eigs, maxEVs, config.problem, config.branch)
    plotSpectrum(eigs, config.plotLabel, config.plotLims)
    # u, v, w = getPerturbation(max_vv, max_eig, var, beta, D1)
    # plotEVector(y_cheb, u, v)
    return

