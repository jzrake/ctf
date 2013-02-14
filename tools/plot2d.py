#!/usr/bin/env python

import sys
import optparse
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14


def plotfile(f, variable, ax):
    h5f = h5py.File(f)
    prim = h5f['prim']
    ax.imshow(prim[variable].value.T, origin='image')


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-o', '--output')
    opts, fnames = parser.parse_args()

    fig = plt.figure(figsize=[10,10])
    ax1 = fig.add_subplot('111')

    for f in fnames:
        plotfile(f, 'rho', ax1)

    if opts.output:
        plt.savefig(opts.output)
    else:
        plt.show()
