import numpy as np
import time
import EVproblem as evp
import postprocessing as pp
import readData as rd
from config import Branch, Config


def runPoiseuille(config: Config) -> complex:
    y_cheb, D, D2, D3, D4 = evp.getChebMatrices(config.n)

    y_cheb_inner, D_inner, D2_inner, D3_inner, D4_inner = evp.clampedBC(
        y_cheb, D, D2, D3, D4
    )


    U = rd.poiseuille_flow(y_cheb_inner)
    U_y = rd.dpoiseuille_flow(y_cheb_inner)
    U_yy = rd.d2poiseuille_flow(y_cheb_inner)

    start = time.time()
    eigs = []
    if config.branch == Branch.Temporal:
        L, M = evp.getOSMatricesTemporal(
            config.n, config.re, config.alpha, config.beta, D2_inner, D4_inner, U, U_yy
        )
        eigs, max_eig, max_vv, _ = evp.solveEVproblemTemporal(L, M, U, U_y, config)
        if config.use_c:
            eigs = eigs / config.alpha
            max_eig = max_eig / config.alpha
    else:
        L, M = evp.getOSMatricesSpatial(
            config.n,
            config.re,
            config.omega,
            config.beta,
            D_inner,
            D2_inner,
            D3_inner,
            D4_inner,
            U,
            U_yy,
        )
        eigs, max_eig, max_vv = evp.solveEVproblemSpatial(L, M, U, U_y, config)
    end = time.time()
    print(f"Time to solve eigenvalue problem: {end - start}")

    if config.doPlot:
        pp.postprocess(eigs, config)

    return max_eig


def runUFromFile(config: Config) -> complex:
    # change of variable parameters
    y_max: float = 25
    yl: float = 5  # position where we want half of the points to be distributed to
    a: float = yl * y_max / (y_max - 2 * yl)
    b: float = 1 + 2 * a / y_max

    start = time.time()
    y_cheb, D, D2, D3, D4 = evp.getChebMatrices(config.n)

    y_cheb_inner, D_inner, D2_inner, D3_inner, D4_inner = evp.clampedBC(
        y_cheb, D, D2, D3, D4
    )

    y_cheb_inner, D_inner, D2_inner, D3_inner, D4_inner = evp.adjustCoordinates(
        y_cheb_inner, D_inner, D2_inner, D3_inner, D4_inner, a, b
    )
    # # U, U_yy = rd.getBlasius(y_cheb_inner, uinf)

    # print("y_cheb_inner", y_cheb_inner)
    # print("D4_inner", D4_inner)

    # print("D_inner", D4_inner)

    # filename = "data/points_x254.dat"

    # rd.blasius_profile(x, filename, uinf, Re)
    y, U, _ = rd.read_baseflow(config.filenameUprofile)

    def df_order_four(f, x):
        h = x[1] - x[0]
        df = np.zeros_like(f)

        # Fourth-order centered difference for interior points
        df[2:-2] = (f[:-4] - 8 * f[1:-3] + 8 * f[3:-1] - f[4:]) / (12 * h)

        # Mixed forward/backward scheme for second and second-to-last points
        df[1] = (-3 * f[0] - 10 * f[1] + 18 * f[2] - 6 * f[3] + f[4]) / (12 * h)
        df[-2] = (3 * f[-1] + 10 * f[-2] - 18 * f[-3] + 6 * f[-4] - f[-5]) / (12 * h)

        # Fourth-order forward difference at the first point
        df[0] = (-25 * f[0] + 48 * f[1] - 36 * f[2] + 16 * f[3] - 3 * f[4]) / (12 * h)

        # Fourth-order backward difference at the last point
        df[-1] = (25 * f[-1] - 48 * f[-2] + 36 * f[-3] - 16 * f[-4] + 3 * f[-5]) / (
            12 * h
        )
        return df

    U_y = df_order_four(U, y)
    U_yy = df_order_four(U_y, y)
    # U_y = np.gradient(U, y, edge_order=2)
    # U_yy = np.gradient(U_y, y, edge_order=2)

    U, U_y, U_yy = evp.interpolateCheb(y_cheb_inner - y_cheb_inner[0], y, U, U_y, U_yy)

    # print("y_cheb_inner", y_cheb_inner)
    # print("U", U)
    # print("U_y", U_y)
    # print("U_yy", U_yy)
    # U_y = D_inner @ U
    # U_yy = D_inner @ U_y

    if config.branch == Branch.Temporal:
        L, M = evp.getOSMatricesTemporal(
            config.n, config.re, config.alpha, config.beta, D2_inner, D4_inner, U, U_yy
        )
        end = time.time()
        print(f"Time to get matrices: {end - start}")
        # print("alpha", alpha)
        # print("beta", beta)

        eigs, max_eig, max_vv, _ = evp.solveEVproblemTemporal(L, M, U, U_y, config)
        start = time.time()
        print(f"Time to solve eigenvalues: {start - end}")

        if config.use_c:
            eigs = eigs / config.alpha
            max_eig = max_eig / config.alpha
    else:
        L, M = evp.getOSMatricesSpatial(
            config.n,
            config.re,
            config.omega,
            config.beta,
            D_inner,
            D2_inner,
            D3_inner,
            D4_inner,
            U,
            U_yy,
        )
        end = time.time()
        print(f"Time to get matrices: {end - start}")

        eigs, max_eig, max_vv = evp.solveEVproblemSpatial(L, M, U, U_y, config)
        start = time.time()
        print(f"Time to solve eigenvalues: {start - end}")

    if config.doPlot:
        pp.postprocess(eigs, config)

    return max_eig

    # def extendBlaisus(filename):
    #     # filename = "data/blasius.dat"
    #     data = np.loadtxt(filename, skiprows=3)
    #     x = data[:, 0]
    #     y = data[:, 1]
    #     U = data[:, 2]
    #     V = data[:, 3]
    #     U_y = data[:, 4]
    #     U_yy = data[:, 5]

    #     y_max = 75
    #     dy = y[1] - y[0]
    #     added_points = np.arange(y[-1] + dy, y_max, dy)
    #     added_x = np.ones(added_points.size) * x[-1]
    #     added_U = np.ones(added_points.size) * U[-1]
    #     added_V = np.ones(added_points.size) * V[-1]
    #     added_U_y = np.ones(added_points.size) * U_y[-1]
    #     added_U_yy = np.ones(added_points.size) * U_yy[-1]

    #     with open(filename, "a") as file:
    #         for i in range(added_points.size):
    #             file.write(
    #                 str(added_x[i])
    #                 + " "
    #                 + str(added_points[i])
    #                 + " "
    #                 + str(added_U[i])
    #                 + " "
    #                 + str(added_V[i])
    #                 + " "
    #                 + str(added_U_y[i])
    #                 + " "
    #                 + str(added_U_yy[i])
    #                 + "\n"
    #             )

    # extendBlaisus(filename)
