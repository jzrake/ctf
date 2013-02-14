#!/usr/bin python

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['xtick.labelsize'] = 14
matplotlib.rcParams['ytick.labelsize'] = 14

L1_Table = {8:{},16:{},32:{},64:{},128:{},256:{},512:{},1024:{}}

def plotfile(f, scheme, ax):
    data = np.loadtxt(f)
    I = np.argsort(data, axis=0)
    N = [ ]
    L = [ ]
    for i in I[:,1]:
        N.append(data[i,0])
        L.append(data[i,1])
    ax.loglog(N, L, '-o', lw=4, ms=10.0, mfc='none', label=scheme)

    for i in range(len(N)-1):
        L1_Table[N[i]][scheme] = -np.log10(L[i+1]/L[i]) / np.log10(N[i+1]/N[i])

if __name__ == "__main__":
    fnames = sys.argv[1:]
    fig = plt.figure()
    ax1 = fig.add_subplot('111')

    for f in fnames:
        scheme = f[:-4].upper()
        plotfile(f, scheme, ax1)

    ax1.xaxis.set_major_formatter(FormatStrFormatter(r'$%d$'))
    ax1.set_xticks([8,16,32,64,128,256,512,1024])
    ax1.set_yticks([1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2])
    ax1.set_xlim(4, 2048)
    ax1.set_xlabel(r'$N$')
    ax1.set_ylabel(r'$L_1$ error')
    ax1.legend(loc='best')
    plt.show()

    schemes = ['HLLC-PLM-MUSCL',
               'HLLC-PLM-RK3',
               'HLLC-WENO5-RK3',
               'CHAR-WENO5-RK3',
               'CHAR-WENO5-RK4']

    print '|$N$|%s|' % '|'.join('=%s='%s for s in schemes)
    for N in sorted(L1_Table.keys())[1:]:
        print '|%d|%s|' % (N, '|'.join(["%3.2f"%L1_Table[N][s] for s in schemes]))

