import matplotlib.pyplot as plt
import os
import numpy as np
from config import Config, Branch, ProblemType
from pathlib import Path


def readData(filename: Path, colX: int, colY: int, numSkipHeaderLines: int):
    data = np.loadtxt(filename, skiprows=numSkipHeaderLines)

    X = data[:, colX]
    Y = data[:, colY]

    return np.array(X), np.array(Y)


def getLabels(conf: Config):
    varlabel = ""
    plotlabel = ""
    if conf.branch == Branch.Temporal:
        varlabel = "α"
        if conf.use_c:
            plotlabel = "c"
        else:
            plotlabel = "ω"
    elif conf.branch == Branch.Spatial:
        varlabel = "ω"
        plotlabel = "α"

    return varlabel, plotlabel


def plotUprofile(conf: Config):
    # Read the data
    y, u = readData(
        conf.filenameUprofile,
        conf.colX,
        conf.colY,
        conf.numSkipHeaderLines,
    )

    filenameBL = "data/blasius.dat"
    pathBL = Path(filenameBL)
    y_bl, u_bl = readData(pathBL, 1, 2, 3)

    u_bl_int = np.interp(y, y_bl, u_bl)

    # print the relative error in Loo norm
    u_diff = np.linalg.norm(u - u_bl_int, ord=np.inf)
    u_diff /= np.linalg.norm(u_bl_int, ord=np.inf)
    print(f"||u_DNS - u_BL||_oo / ||u_BL||_oo = {u_diff:.2e}")

    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(u, y, label="Blaisus DNS", alpha=0.7)
    ax.plot(u_bl, y_bl, label="Blasius", alpha=0.7)
    ax.set_ylabel("y")
    ax.set_xlabel("U")
    ax.legend()

    plt.show()


def plotEV(conf: Config):
    re, im = readData(conf.fileWriteEigenvalues, 0, 1, 0)
    varlabel, plotlabel = getLabels(conf)
    if conf.problem == ProblemType.BoundaryLayer:
        conf.vars_range_r *= conf.DELTASTAR
        conf.vars_range_i *= conf.DELTASTAR

    if conf.run_multiple:
        # plot all the points with same imaginary part with same color
        # check how many different imaginary parts are there

        # vars_r = np.repeat(conf.vars_range_r, len(conf.vars_range_i))
        vars_i = np.tile(conf.vars_range_i, len(conf.vars_range_r))

        if len(conf.vars_range_i) > 1:
            for im_val in conf.vars_range_i:
                mask = vars_i == im_val
                plt.scatter(
                    re[mask],
                    im[mask],
                    label=f"Im({varlabel}) = {np.round(im_val, 3)}",
                    alpha=0.7,
                )
        else:
            # rainbow colors
            sc = plt.scatter(
                re,
                im,
                c=conf.vars_range_r,
                cmap="rainbow",
                label=f"Im({varlabel}) = {np.round(conf.vars_range_i[0], 3)}",
            )

            plt.colorbar(sc, label=f"Re({varlabel})")
        plt.xlim(conf.plot_lims.xmin, conf.plot_lims.xmax)
        plt.ylim(conf.plot_lims.ymin, conf.plot_lims.ymax)
    else:
        plt.scatter(re, im, label="Eigenvalues", alpha=0.7)
        plt.xlim(conf.plot_lims.xmin, conf.plot_lims.xmax)
        plt.ylim(conf.plot_lims.ymin, conf.plot_lims.ymax)

    plt.xlabel(f"Re({plotlabel})")
    plt.ylabel(f"Im({plotlabel})")
    plt.title("Eigenvalues")
    plt.grid()
    plt.legend()
    plt.show()

    return


def main():
    # read data
    CONFIG_FILE = "config/input.toml"
    DIR_SCRIPT = os.path.dirname(os.path.realpath(__file__))
    CONFIG_FILE = os.path.join(DIR_SCRIPT, "../", CONFIG_FILE)

    conf = Config.from_toml(CONFIG_FILE)

    if conf.problem == ProblemType.Custom and conf.plotUprofile:
        plotUprofile(conf)

    plotEV(conf)


if __name__ == "__main__":
    main()
