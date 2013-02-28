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

exact = True
problem = ""
outtime = None

def plotfile(f, variable, ax):
    global problem, outtime
    h5f = h5py.File(f)
    problem = h5f["problem"].value
    outtime = h5f["time"].value
    prim = h5f['prim' ].value
    exac = h5f['exact'].value
    ivar = {'rho': 0, 'pre': 1, 'vx': 2, 'vy': 3, 'vz': 4}[variable]
    lvar = {'rho': 'density', 'pre': 'pressure', 'vx': 'velocity'}[variable]
    x = np.linspace(0, 1, prim.shape[0])
    if exact:
        ax.plot(x, exac[:,ivar],  'o', mfc='none', label='EXACT')
    ax.plot(x, prim[:,ivar], mfc='none', lw=1.5, label=h5f['id'].value.upper())
    ax.set_ylabel(lvar)


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-o', '--output')
    opts, fnames = parser.parse_args()

    fig = plt.figure(figsize=[8,10])
    ax1 = fig.add_subplot('311')
    ax2 = fig.add_subplot('312')
    ax3 = fig.add_subplot('313')

    for f in fnames:
        plotfile(f, 'rho', ax1)
        plotfile(f, 'pre', ax2)
        plotfile(f, 'vx' , ax3)
        exact = False

    for ax in [ax1,ax2,ax3]:
        yl = list(ax.get_ylim())
        yl[0] += 1e-6
        yl[1] -= 1e-6
        ax.set_ylim(yl)

    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xlabel(r'position $x$')

    fig.suptitle(r"$\rm{%s} \qquad t=%3.2f$" % (problem, outtime), fontsize=18)
    plt.subplots_adjust(hspace=0)
    ax1.legend(loc='best')

    if opts.output:
        plt.savefig(opts.output)
    else:
        plt.show()
