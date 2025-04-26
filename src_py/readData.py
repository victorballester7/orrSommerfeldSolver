import numpy as np
import numpy.typing as npt
from typing import Tuple
from pathlib import Path
import sys
import os

# Add boeingGap to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from src.blasiusIncNS import solve_bvp as sbvp
from src.blasiusIncNS import physicalQuantities as pq

def read_baseflow(filename: Path):
    data = np.loadtxt(filename, skiprows=3)
    y = data[:, 1]
    U = data[:, 2]
    V = data[:, 3]
    # U_y = data[:, 4]
    # U_yy = data[:, 5]
    return y, U, V

def poiseuille_flow(y: npt.NDArray[np.float64]):
    return 1 - y**2

def dpoiseuille_flow(y: npt.NDArray[np.float64]):
    return -2 * y

def d2poiseuille_flow(y: npt.NDArray[np.float64]):
    return -2 * np.ones(y.size)


def writeSol2File(x: npt.NDArray[np.float64], y: npt.NDArray[np.float64], filename: Path, u_inf: float, re_x: float):
    with open(filename, "w") as file:
        file.write("# Blasius solution\n")
        file.write("\n")
        file.write("# x y U V\n")
        for i in range(x.size):
            u = u_inf * y[1][i]
            u_y = u_inf * y[2][i]
            u_yy = u_inf * (- 0.5 * y[0][i] * y[2][i])
            v = 0.5 * u_inf / np.sqrt(re_x) * (x[i] * y[1][i] - y[0][i])
            file.write("0 " + str(x[i]) + " " + str(u) + " " + str(v) + " " + str(u_y) + " " + str(u_yy) + "\n")

def blasius_profile(x: npt.NDArray[np.float64], filename: Path, u_inf: float, re_deltaStar: float):
    solution = sbvp.solve_BVP(3, x)
    y = solution.sol(x)
    re_x = pq.getRe_x(re_deltaStar,  pq.computeDeltaStar(y[0], x[-1]))
    writeSol2File(x, y, filename, u_inf, re_x)
    return

def getBlasius(x: npt.NDArray[np.float64], u_inf: float) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    solution = sbvp.solve_BVP(3, x)
    y = solution.sol(x)
    u = u_inf * y[1]
    u_yy = u_inf * (- 0.5 * y[0] * y[2])
    return u, u_yy

# def addPoints(y, U, V, max_y):
def addPoints(y: npt.NDArray[np.float64], U: npt.NDArray[np.float64], V: npt.NDArray[np.float64], max_y: float) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    dy = y[-1] - y[-2]
    added_points = np.arange(y[-1] + dy, max_y, dy)
    added_U = np.ones(added_points.size) * U[-1]
    added_V = np.ones(added_points.size) * V[-1]

    y = np.append(y, added_points)
    U = np.append(U, added_U)
    V = np.append(V, added_V)
    return y, U, V


def readEValues(filename: str):
    data = np.loadtxt(filename)
    reals = data[:, 0]
    imags = data[:, 1]
    e_values = reals + 1j * imags
    return e_values
