import ROOT
import numpy as np
import h5py
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

f=h5py.File('MagFieldData_array_5cm_grid_gpu_v3_debug.hdf5','r')
array_bx = f['Bx_map']
array_by = f['By_map']
array_bz = f['Bz_map']

plot=plt.figure()
ax=plot.gca()

x = np.linspace(-500.0,500.0, num=201)
y = np.linspace(-500.0,500.0, num=201)
z = np.linspace(-500.0,500.0, num=201)

Bx = []
By = []
Bz = []

#print(array_bx[1.0,1.0,1.0])
#print(array_by[1.0,1.0,1.0])
#print(array_bz[40.0,40.0,50.0])

for xi in range(-500,501,5):
	bx=[]
	by=[]
	bz=[]
	for zk in range(-500,501,5):
		#print(str(array_bx[xi,0,zk])+'     '+str(array_bz[xi,0,zk]))
		bx.append(array_bx[xi,0,zk])
		by.append(array_by[xi,0,zk])
		bz.append(array_bz[xi,0,zk])
		B = np.sqrt(array_bx[xi,0,zk]*array_bx[xi,0,zk] + array_bz[xi,0,zk]*array_bz[xi,0,zk])
		if (B>15.0):
			print(str(xi)+'     '+str(zk)+'     '+str(array_bx[xi,0,zk])+'      '+str(array_by[xi,0,zk])+'      '+str(array_bz[xi,0,zk]))
	Bx.append(bx)
	By.append(by)
	Bz.append(bz)

_Bx=np.asarray(Bx)
_By=np.asarray(By)
_Bz=np.asarray(Bz)

B = np.sqrt(_Bx*_Bx + _By*_By+ _Bz*_Bz)
lw= 5.0*_Bz/B.max()
print(B.max())

#if B == B.max():
	
strm=ax.streamplot(x,z,_Bx,_Bz, density=[0.5,1.0], color=B, linewidth=lw)
plot.colorbar(strm.lines)
plt.show(plot)
plt.savefig('stream.png')
