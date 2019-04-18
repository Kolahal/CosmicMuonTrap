import numpy as np
import h5py
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

f=h5py.File('MagFieldData_list_10m.hdf5','r')
Ls = f['B_map']

soa = np.array(Ls)
X, Y, Z, U, V, W = zip(*soa)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(X, Y, Z, U, V, W, length=0.01)
ax.set_xlim([-50, 50])
ax.set_ylim([-50, 50])
ax.set_zlim([0., 100])
plt.show()
plt.savefig('abc.png')

