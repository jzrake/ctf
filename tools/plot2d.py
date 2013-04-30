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
    ax.imshow(prim[variable].value.T, origin='image', interpolation='nearest')
    return prim[variable].shape[0]


def plotfile2(f, variable, ax):
    h5f = h5py.File(f)
    prim = h5f['prim']
    vx = prim['vx'].value.T
    vy = prim['vy'].value.T
    x = np.linspace(-0.5, 0.5, vx.shape[0])
    y = np.linspace(-0.5, 0.5, vy.shape[1])
    ax.streamplot(x, y, vx, vy)
    return prim[variable].shape[0]


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-o', '--output')
    opts, fnames = parser.parse_args()

    fig = plt.figure(figsize=[10,10])
    ax1 = fig.add_subplot('221')
    ax2 = fig.add_subplot('222')
    ax3 = fig.add_subplot('223')
    ax4 = fig.add_subplot('224')
    axs = [ax1, ax2, ax3, ax4]

    for ax,f in zip(axs,fnames):
        N = plotfile(f, 'pre', ax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.02, 0.02, r'$%d^2$'%N, fontsize=24,
                color='white',
                transform=ax.transAxes,)

    plt.subplots_adjust(left=0.03, bottom=0.03, right=0.98, top=0.98,
                        wspace=0, hspace=0)
    if opts.output:
        plt.savefig(opts.output)
    else:
        plt.show()
