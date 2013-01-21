import h5py
import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot('411')
ax2 = fig.add_subplot('412')
ax3 = fig.add_subplot('413')
ax4 = fig.add_subplot('414')

h5f = h5py.File("euler.h5")
ax1.plot(h5f["prim"][:,0])
ax2.plot(h5f["prim"][:,1])
ax3.plot(h5f["grav"][:,0])
ax4.plot(h5f["grav"][:,1])

plt.show()
