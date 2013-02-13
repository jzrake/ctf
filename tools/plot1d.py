
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt


def plotfile(f, variable, ax):
    h5f = h5py.File(f)
    prim = h5f['prim' ].value
    exac = h5f['exact'].value
    ivar = {'rho': 0, 'pre': 1, 'vx': 2, 'vy': 3, 'vz': 4}[variable]
    ax.plot(prim[:,ivar])
    ax.plot(exac[:,ivar])
    ax.set_ylabel(variable)


if __name__ == "__main__":
    fnames = sys.argv[1:]

    fig = plt.figure()
    ax1 = fig.add_subplot('311')
    ax2 = fig.add_subplot('312')
    ax3 = fig.add_subplot('313')

    for f in fnames:
        plotfile(f, 'rho', ax1)
        plotfile(f, 'pre', ax2)
        plotfile(f, 'vx' , ax3)

    plt.show()
