
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt

exact = True
problem = ""

def plotfile(f, variable, ax):
    global problem
    h5f = h5py.File(f)
    problem = h5f["problem"].value
    prim = h5f['prim' ].value
    exac = h5f['exact'].value
    ivar = {'rho': 0, 'pre': 1, 'vx': 2, 'vy': 3, 'vz': 4}[variable]
    lvar = {'rho': 'density', 'pre': 'pressure', 'vx': 'velocity'}[variable]
    if exact:
        ax.plot(exac[:,ivar],  'o', mfc='none', label='EXACT')
    ax.plot(prim[:,ivar], mfc='none', lw=1.5, label=h5f['id'].value.upper())
    ax.set_ylabel(lvar)
    ax.set_xlim([-2, prim.shape[0]+2])

if __name__ == "__main__":
    fnames = sys.argv[1:]

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
        ax.set_ylim([ax.get_ylim()[0]/1.1, ax.get_ylim()[1]*1.1])

    ax1.set_xticks([])
    ax2.set_xticks([])

    fig.suptitle(problem)
    plt.subplots_adjust(hspace=0)
    ax1.legend(loc='best')
    plt.show()
