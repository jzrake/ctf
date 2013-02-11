
import glob
import numpy as np
import h5py
import matplotlib.pyplot as plt
import streamplot

def show_frame(fname, show=True):
    h5f = h5py.File(fname)
    vx = h5f['prim'][...,2]
    vy = h5f['prim'][...,3]
    Bx = h5f['prim'][...,5]
    By = h5f['prim'][...,6]
    X = np.linspace(-1, 1, Bx.shape[0])
    Y = np.linspace(-1, 1, Bx.shape[1])
    kwargs = dict(arrowsize=1,
                  density=1,
                  linewidth=1)
    streamplot.streamplot(X, Y, vx.T, vy.T, color='b', **kwargs)
    streamplot.streamplot(X, Y, Bx.T, By.T, color='r', **kwargs)
    plt.axis('equal')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    if show:
        plt.show()
    else:
        plt.savefig(fname.replace('.h5', '.png'))
        plt.clf()


def show_density(fname, show=True):
    h5f = h5py.File(fname)
    rho = h5f['prim'][...,0].T
    plt.imshow(rho, origin='image')
    plt.colorbar()
    plt.axis('equal')
    if show:
        plt.show()
    else:
        plt.savefig(fname.replace('.h5', '.png'))
        plt.clf()


fnames = glob.glob('data/*.h5')
show_density(fnames[-1])
#for fname in fnames:
#    print fname
#    show_frame(fname, show=False)
