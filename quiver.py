import numpy as np
import h5py
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

'''
f=h5py.File('MagFieldData_list_3D_gpu.hdf5','r')
Ls = f['B_map']

soa = np.array(Ls)
X, Y, Z, U, V, W = zip(*soa)

fig = plt.figure()
'''

soa = np.array([[0, 0, 3, 2], [0, 0, 1, 1], [0, 0, 9, 9]])
X, Y, U, V = zip(*soa)

plt.savefig('abc.png')
plt.figure()
ax = plt.gca()
ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1)
ax.set_xlim([-1, 10])
ax.set_ylim([-1, 10])
plt.draw()
plt.show()
plt.savefig('abc.png')

'''
ax = fig.add_subplot(111, projection='3d')
ax.quiver(X, Y, Z, U, V, W, length=0.01)
ax.set_xlim([-50, 50])
ax.set_ylim([-50, 50])
ax.set_zlim([0., 100])
plt.show()
plt.savefig('abc.png')
'''
