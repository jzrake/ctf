import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt

def gravity():
    fig = plt.figure()
    ax1 = fig.add_subplot('411')
    ax2 = fig.add_subplot('412')
    ax3 = fig.add_subplot('413')
    ax4 = fig.add_subplot('414')

    dx = 1.0 / 94.0
    h5f = h5py.File("euler.h5")
    rho = h5f["prim"][:,0]
    phi = h5f["grav"][:,0]
    gph = (np.roll(phi, -1) - np.roll(phi, +1)) / (2*dx)
    lph = (np.roll(gph, -1) - np.roll(gph, +1)) / (2*dx)

    ax1.plot(h5f["prim"][:,0], label='rho')
    ax1.plot(lph+rho.mean(), label='laplacian phi + rhobar')
    ax1.legend()
    ax2.plot(h5f["prim"][:,1])
    ax3.plot(h5f["grav"][:,0], label='')
    ax4.plot(h5f["grav"][:,1], label='spectral')
    ax4.plot(gph, label='differenced')
    ax4.legend()
    plt.show()

def hydro():
    fig = plt.figure()
    ax1 = fig.add_subplot('311')
    ax2 = fig.add_subplot('312')
    ax3 = fig.add_subplot('313')

    for fname in glob.glob('data/*.h5')[-1:]:
        h5f = h5py.File(fname)
        ax1.plot(h5f["prim" ][:,0], label='code')
        ax1.plot(h5f["exact"][:,0], label='exact')
        ax2.plot(h5f["prim" ][:,1], label='code')
        ax2.plot(h5f["exact"][:,1], label='exact')
        ax3.plot(h5f["prim" ][:,2], label='code')
        ax3.plot(h5f["exact"][:,2], label='exact')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    plt.show()

hydro()
