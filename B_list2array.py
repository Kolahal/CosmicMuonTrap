import numpy as np
import h5py
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import time

f=h5py.File('MagFieldData_list_3D.hdf5','r')
Ls = f['B_map']

print(str(len(Ls)))
Arr_bx = np.zeros((101,101,101))
Arr_by = np.zeros((101,101,101))
Arr_bz = np.zeros((101,101,101))
tic = time.time()
for i in range(len(Ls)):
	Arr_bx[int(Ls[i][0]),int(Ls[i][1]),int(Ls[i][2])] = Ls[i][3]
	Arr_by[int(Ls[i][0]),int(Ls[i][1]),int(Ls[i][2])] = Ls[i][4]
	Arr_bz[int(Ls[i][0]),int(Ls[i][1]),int(Ls[i][2])] = Ls[i][5]

with h5py.File('MagFieldData_array_3D.hdf5','w') as hf:
	hf.create_dataset('Bx_map',data=Arr_bx)
	hf.create_dataset('By_map',data=Arr_by)
	hf.create_dataset('Bz_map',data=Arr_bz)

toc = time.time()
print(str(toc-tic))
