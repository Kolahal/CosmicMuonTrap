import numpy as np
import h5py
from numpy import genfromtxt
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.gridspec as gridspec
import time

f=h5py.File('MagFieldData_list_3D_gpu_10m.hdf5','r')
Ls = f['B_map']

print(str(len(Ls)))
Arr_bx = np.zeros((501,501,501))
Arr_by = np.zeros((501,501,501))
Arr_bz = np.zeros((501,501,501))

print(Arr_bz.shape)

tic = time.time()
for i in range(len(Ls)):
	#print(str(i)+'    '+str(Ls[i][0])+'    '+str(Ls[i][1])+'   '+str(Ls[i][2]))
	Arr_bx[int(Ls[i][0]),int(Ls[i][1]),int(Ls[i][2])] = Ls[i][3]
	Arr_by[int(Ls[i][0]),int(Ls[i][1]),int(Ls[i][2])] = Ls[i][4]
	Arr_bz[int(Ls[i][0]),int(Ls[i][1]),int(Ls[i][2])] = Ls[i][5]
	#print(str(i)+'    '+str(Ls[i][3])+'    '+str(Ls[i][4])+'   '+str(Ls[i][5]))

Arr_bx = Arr_bx[np.nonzero(Arr_bx)]
Arr_by = Arr_by[np.nonzero(Arr_by)]
Arr_bz = Arr_bz[np.nonzero(Arr_bz)]

print(Arr_bz.shape)

with h5py.File('MagFieldData_array_3D_gpu_10m.hdf5','w') as hf:
	hf.create_dataset('Bx_map',data=Arr_bx)
	hf.create_dataset('By_map',data=Arr_by)
	hf.create_dataset('Bz_map',data=Arr_bz)

toc = time.time()
print(str(toc-tic))
