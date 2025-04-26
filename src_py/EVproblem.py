import numpy as np
import numpy.typing as npt
import scipy.linalg as scla
import dmsuite as dm
import postprocessing as pp
from config import Config
from typing import Tuple


def mapFromCheb(
    y: npt.NDArray[np.float64], a: float, b: float
) -> npt.NDArray[np.float64]:
    # mapping taken from:
    # Numerical Methods for Hypersonic Boundary-Layer Stability
    # M. R. MALIK
    return a * (1.0 + y) / (b - y)


def mapToCheb(
    y_cheb: npt.NDArray[np.float64], a: float, b: float
) -> npt.NDArray[np.float64]:
    return (y_cheb * b - a) / (y_cheb + a)


def dmapToCheb(
    y_cheb: npt.NDArray[np.float64], a: float, b: float
) -> npt.NDArray[np.float64]:
    diff = a * (b + 1) / (a + y_cheb) ** 2
    return diff.astype(np.float64)


def d2mapToCheb(
    y_cheb: npt.NDArray[np.float64], a: float, b: float
) -> npt.NDArray[np.float64]:
    diff = -2.0 * a * (b + 1) / (a + y_cheb) ** 3
    return diff.astype(np.float64)


def d3mapToCheb(
    y_cheb: npt.NDArray[np.float64], a: float, b: float
) -> npt.NDArray[np.float64]:
    diff = 6.0 * a * (b + 1) / (a + y_cheb) ** 4
    return diff.astype(np.float64)


def d4mapToCheb(
    y_cheb: npt.NDArray[np.float64], a: float, b: float
) -> npt.NDArray[np.float64]:
    diff = -24.0 * a * (b + 1) / (a + y_cheb) ** 5
    return diff.astype(np.float64)


def adjustCoordinates(
    y_cheb: npt.NDArray[np.float64],
    D: npt.NDArray[np.float64],
    D2: npt.NDArray[np.float64],
    D3: npt.NDArray[np.float64],
    D4: npt.NDArray[np.float64],
    a: float,
    b: float,
) -> Tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]:
    # adjust coordinates of the derivatives due to change of variable formula
    y_cheb_new = mapFromCheb(y_cheb, a, b)
    dmap = dmapToCheb(y_cheb_new, a, b)
    d2map = d2mapToCheb(y_cheb_new, a, b)
    d3map = d3mapToCheb(y_cheb_new, a, b)
    d4map = d4mapToCheb(y_cheb_new, a, b)

    D_new = np.diag(dmap) @ D
    D2_new = np.diag(dmap**2) @ D2 + np.diag(d2map) @ D
    D3_new = (
        np.diag(dmap**3) @ D3 + 3.0 * np.diag(dmap * d2map) @ D2 + np.diag(d3map) @ D
    )
    D4_new = (
        np.diag(dmap**4) @ D4
        + 6.0 * np.diag(dmap**2 * d2map) @ D3
        + np.diag(3.0 * d2map**2 + 4 * dmap * d3map) @ D2
        + np.diag(d4map) @ D
    )

    return y_cheb_new, D_new, D2_new, D3_new, D4_new


def interpolateCheb(
    y_cheb: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    U: npt.NDArray[np.float64],
    U_y: npt.NDArray[np.float64],
    U_yy: npt.NDArray[np.float64],
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    # a = y[0]
    # b = y[-1]
    # y_toCheb = Real2Cheb(a, b, y)
    # U_toCheb = Real2Cheb(a, b, U)
    # U_yy_toCheb = Real2Cheb(a, b, U_yy)
    U_new = np.interp(y_cheb, y, U)
    U_y_new = np.interp(y_cheb, y, U_y)
    U_yy_new = np.interp(y_cheb, y, U_yy)
    return U_new, U_y_new, U_yy_new


def getChebMatrices(
    N: int,
) -> Tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    y_cheb, D = dm.chebdif(N, 1)
    D = D[0, :, :]

    # flip back the vector and the matrix
    y_cheb = y_cheb[::-1]
    D = D[::-1, ::-1]

    D2 = D @ D
    D3 = D2 @ D
    D4 = D2 @ D2
    return y_cheb, D, D2, D3, D4


def clampedBC(
    y_cheb: npt.NDArray[np.float64],
    D1: npt.NDArray[np.float64],
    D2: npt.NDArray[np.float64],
    D3: npt.NDArray[np.float64],
    D4: npt.NDArray[np.float64],
) -> Tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]:
    y_cheb2 = y_cheb * y_cheb
    s = np.hstack([0.0, 1.0 / (1.0 - y_cheb2[1:-1]), 0.0])
    # print("s")
    # print(np.diag(1.0 - y_cheb2) @ D4 + np.diag(-8.0 * y_cheb) @ D3 - 12.0 * D2)

    D4_new = (
        np.diag(1.0 - y_cheb2) @ D4 + np.diag(-8.0 * y_cheb) @ D3 - 12.0 * D2
    ) @ np.diag(s)
    D1 = D1[1:-1, 1:-1]
    D2 = D2[1:-1, 1:-1]
    D3 = D3[1:-1, 1:-1]
    D4 = D4_new[1:-1, 1:-1]
    return y_cheb[1:-1], D1, D2, D3, D4


def getOSMatricesTemporal(
    N: int,
    Re: float,
    alpha: complex,
    beta: complex,
    D2: npt.NDArray[np.float64],
    D4: npt.NDArray[np.float64],
    U: npt.NDArray[np.float64],
    U_yy: npt.NDArray[np.float64],
) -> Tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    id = np.identity(N - 2)

    # Orr-Sommerfeld eigenvalue problem:
    # L * v = omega * (M * v)

    # Squire's eq.
    # A * eta = B * v

    k2 = alpha**2 + beta**2

    # The matrices are those on the interior points of the domain. The BCs are alrewady implemented in D4 (neumann) but not in D2 because of:
    # (
    # taken from the paper of
    # A MATLAB Differentiation Matrix Suite
    # J. A. C. WEIDEMAN
    # University of Stellenbosch
    # and
    # S. C. REDDY
    # Oregon State University
    # )
    # Let D4 be the fourth-derivative Chebyshev matrix that implements the clamped (i.e. y(\pm 1) = y'(\pm 1) = 0) boundary conditions, and let D2 be the second-derivative Chebyshev matrix with boundary conditions y(\pm 1) = 0. Note that different boundary conditions are employed for these matrices—this approach eliminates spurious eigenvalues [Huang and Sloan 1994].
    # the dirichlet BCs are naturally impleemented on both D2 and D4 by considering only integration on the interior points of the domain

    L = (D4 - 2 * k2 * D2 + k2**2 * id) / Re + 1j * alpha * (
        np.diag(U_yy) - np.diag(U) @ (D2 - id * k2)
    )
    M = 1j * (id * k2 - D2)

    return L, M


def getOSMatricesSpatial(
    N: int,
    Re: float,
    omega: complex,
    beta: complex,
    D: npt.NDArray[np.float64],
    D2: npt.NDArray[np.float64],
    D3: npt.NDArray[np.float64],
    D4: npt.NDArray[np.float64],
    U: npt.NDArray[np.float64],
    U_yy: npt.NDArray[np.float64],
) -> Tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    # check page 259 of Schmid and Henningson
    id = np.identity(N - 2)
    zero = np.zeros((N - 2, N - 2))

    # Orr-Sommerfeld eigenvalue problem:
    # L * [alpha * v, v, eta]^T = alpha * M * [alpha * v, v, eta]^T
    # L = [ -R1  -R0  0 ]
    #     [  I    0   0 ]
    #     [  0   -S  -T0]
    # M = [  R2   0   0 ]
    #     [  0    I   0 ]
    #     [  0    0   T1]

    R0 = (
        1j * omega * D2
        + 1.0 / Re * D4
        - 1j * omega * beta**2 * id
        - 2.0 / Re * beta**2 * D2
        + 1.0 / Re * beta**4 * id
    )
    R1 = (
        -2j * omega * D
        - 4.0 / Re * D3
        - 1j * np.diag(U) @ D2
        + 1j * np.diag(U_yy)
        + 4.0 / Re * beta**2 * D
        + 1j * beta**2 * np.diag(U) @ id
    )
    R2 = 4.0 / Re * D2 + 2j * np.diag(U) @ D
    # T1 = 2.0 / Re * D + 1j * np.diag(U)
    # T0 = -1j * omega * id - 1.0 / Re * D2 + 1j * beta**2 / Re * id
    # S = 1j * beta * np.diag(U_y)

    L = np.block(
        [
            [-R1, -R0],
            [id, zero],
        ]
    )
    M = np.block(
        [
            [R2, zero],
            [zero, id],
        ]
    )

    return L, M


def solveEVproblemTemporal(
    L: npt.NDArray[np.complex128],
    M: npt.NDArray[np.complex128],
    U: npt.NDArray[np.float64],
    U_y: npt.NDArray[np.float64],
    config: Config,
) -> Tuple[
    npt.NDArray[np.complex128],
    complex,
    npt.NDArray[np.complex128],
    npt.NDArray[np.complex128],
]:
    # Orr-Sommerfeld eq.
    M_inv = scla.inv(M)
    L = M_inv @ L

    solution = scla.eig(L)

    omega = solution[0]
    vv = solution[1].astype(np.complex128)

    max_omega, max_vv = pp.getMostUnstableEV(omega, vv, config.problem, config.branch)
    # Squire's eq.

    dim = L.shape[0]

    A = max_omega * np.identity(dim) - config.alpha * np.diag(U) - M / config.re
    B = config.beta * np.diag(U_y)

    A_inv = scla.inv(A)
    max_eta = A_inv @ B @ max_vv

    # this method is slower and less accurate
    # omega, vv = scla.eig(L, M)
    return omega, max_omega, max_vv, max_eta


def solveEVproblemSpatial(
    L: npt.NDArray[np.complex128],
    M: npt.NDArray[np.complex128],
    U: npt.NDArray[np.float64],
    U_y: npt.NDArray[np.float64],
    config: Config,
) -> Tuple[npt.NDArray[np.complex128], complex, npt.NDArray[np.complex128]]:
    # Orr-Sommerfeld eq.
    M_inv = scla.inv(M)
    L = M_inv @ L
    solution = scla.eig(L)
    alphas = solution[0]
    vv = solution[1].astype(np.complex128)

    dim = int(L.shape[0])
    vv = vv[dim // 2 :]

    max_alpha, max_vv = pp.getMostUnstableEV(alphas, vv, config.problem, config.branch)
    # Squire's eq.

    # A = max_omega * np.identity(dim) - alpha * np.diag(U) - M / Re
    # B = beta * np.diag(U_y)

    # A_inv = scla.inv(A)
    # max_eta = A_inv @ B @ max_vv

    # this method is slower and less accurate
    # omega, vv = scla.eig(L, M)
    return alphas, max_alpha, max_vv
