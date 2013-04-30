#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py

try:
    h5f = h5py.File(sys.argv[1], 'r')
except:
    print "usage: streamplot infile.h5"
    exit()

def draw_streamlines(fx, fy, title=None):
    plt.figure()
    plt.axis('equal')
    plt.title(title)
    nz = h5f['prim/vx'].shape[2]
    vx = h5f['prim/%s'%fx][:,:,nz/2].T
    vy = h5f['prim/%s'%fy][:,:,nz/2].T
    scalar = vx**2 + vy**2
    x = np.linspace(-0.5, 0.5, vx.shape[0])
    y = np.linspace(-0.5, 0.5, vx.shape[1])
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5, 0.5)
    plt.streamplot(x, y, vx, vy, color=scalar, density=4, cmap='cool')


def draw_slice():
    plt.figure()
    ny = h5f['prim/vz'].shape[1]
    vz = h5f['prim/vz'][:,ny/2,:].T
    plt.imshow(vz, interpolation='nearest')


def draw_line_x():
    plt.figure()
    nx, ny, nz = h5f['prim/vx'].shape
    vx = h5f['prim/vx'][nx/2,:,nz/2]
    vy = h5f['prim/vy'][nx/2,:,nz/2]
    vz = h5f['prim/vz'][nx/2,:,nz/2]
    P = h5f['prim/pre'][nx/2,:,nz/2]
    plt.plot(vx, label='vx')
    plt.plot(vy, label='vy')
    plt.plot(vz, label='vz')
    plt.plot(P, label='P')
    plt.legend()


def draw_line_z():
    plt.figure()
    nx, ny, nz = h5f['prim/vx'].shape
    vx = h5f['prim/vx'][nx/2,ny/2,:]
    vy = h5f['prim/vy'][nx/2,ny/2,:]
    vz = h5f['prim/vz'][nx/2,ny/2,:]
    P = h5f['prim/pre'][nx/2,ny/2,:]
    plt.plot(vx, label='vx')
    plt.plot(vy, label='vy')
    plt.plot(vz, label='vz')
    plt.plot(P, label='P')
    plt.legend()



draw_streamlines('vx', 'vy', title="velocity")
draw_streamlines('Bx', 'By', title="magnetic")
draw_line_z()
draw_slice()
plt.show()
