import sys
import scipy as sp
import scipy.special
import numpy as np
import h5py
import itertools
from itertools import product
import time
from math import pi
import multiprocessing
from multiprocessing import Pool #  Process pool
from multiprocessing import sharedctypes


MMT_dim = int(sys.argv[1])     # m
mu_0    = 4*pi*1.e-7            # Tm/A
a       = []                    # m
I       = []                    # A
C       = []                    # (Tm/A * A)

for zi in range(0,10*MMT_dim+1):
        a.append(0.05)
        I.append(1.e5)
        #print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
        C.append(mu_0*(I[zi]/pi))
for zi in range(10*MMT_dim+1,25*MMT_dim+1):
        a.append(0.75)
        I.append(5.e4)
        #print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
        C.append(mu_0*(I[zi]/pi))
for zi in range(25*MMT_dim+1,75*MMT_dim+1):
        a.append(1.00)
        I.append(1.e4)
        #print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
        C.append(mu_0*(I[zi]/pi))
for zi in range(75*MMT_dim+1,90*MMT_dim+1):
        a.append(0.75)
        I.append(5.e4)
        #print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
        C.append(mu_0*(I[zi]/pi))
for zi in range(90*MMT_dim+1,100*MMT_dim+1):
        a.append(0.05)
        I.append(1.e5)
        #print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
        C.append(mu_0*(I[zi]/pi))

XX = range(-50*MMT_dim,50*MMT_dim+1)
YY = range(-50*MMT_dim,50*MMT_dim+1)
ZZ = range(100*MMT_dim+1)

tic = time.time()

paramlist = list(itertools.product(XX,YY,ZZ))

result_bx = np.ctypeslib.as_ctypes(np.zeros((101,101,101)))
shared_array_bx = sharedctypes.RawArray(result_bx._type_, result_bx)

result_by = np.ctypeslib.as_ctypes(np.zeros((101,101,101)))
shared_array_by = sharedctypes.RawArray(result_by._type_, result_by)

result_bz = np.ctypeslib.as_ctypes(np.zeros((101,101,101)))
shared_array_bz = sharedctypes.RawArray(result_bz._type_, result_bz)

def calculateField(paramlist):
	xi= paramlist[0]
	yj= paramlist[1]
	zk= paramlist[2]
	
	tmp_bx = np.ctypeslib.as_array(shared_array_bx)
	tmp_by = np.ctypeslib.as_array(shared_array_by)
	tmp_bz = np.ctypeslib.as_array(shared_array_bz)
	
	x = 1.e-2*xi
	y = 1.e-2*yj
	
	_bx = 0.0
	_by = 0.0
	_bz = 0.0
	
	for z_n in ZZ:          #iterate over the current loops
		z = 1.e-2*(zk - z_n)
		
		if x==0.0:
			x = 1.e-3
		if y==0.0:
			y = 1.e-3
		if z==0.0:
			z = 1.e-3
		
		r       = np.sqrt(x*x + y*y + z*z)
		rho     = np.sqrt(x*x + y*y)
		alpha   = np.sqrt(a[z_n]*a[z_n] + r*r - 2*a[z_n]*rho)
		beta    = np.sqrt(a[z_n]*a[z_n] + r*r + 2*a[z_n]*rho)
		k       = np.sqrt(1.0-((alpha*alpha)/(beta*beta)))
		
		_bx += ((C[z_n]*x*z)/(2*alpha*alpha*beta*rho*rho))*((a[z_n]*a[z_n]+r*r)*scipy.special.ellipe(k*k) - alpha*alpha*scipy.special.ellipk(k*k))
		_by += ((C[z_n]*y*z)/(2*alpha*alpha*beta*rho*rho))*((a[z_n]*a[z_n]+r*r)*scipy.special.ellipe(k*k) - alpha*alpha*scipy.special.ellipk(k*k))
		_bz += (C[z_n]/(2*alpha*alpha*beta))*((a[z_n]*a[z_n]-r*r)*scipy.special.ellipe(k*k) + alpha*alpha*scipy.special.ellipk(k*k))
		
	tmp_bx[xi, yj, zk] = _bx
	tmp_by[xi, yj, zk] = _by
	tmp_bz[xi, yj, zk] = _bz
	#print(str(xi)+'     '+str(yj)+'     '+str(zk)+'     '+str(_bx)+'     '+str(_by)+'     '+str(_bz))

p = Pool()
res = p.map(calculateField, paramlist)
result_bx = np.ctypeslib.as_array(shared_array_bx)
result_by = np.ctypeslib.as_array(shared_array_by)
result_bz = np.ctypeslib.as_array(shared_array_bz)

with h5py.File('MagFieldData_array_3D.hdf5','w') as hf:
	hf.create_dataset('Bx_map',data=result_bx)
	hf.create_dataset('By_map',data=result_by)
	hf.create_dataset('Bz_map',data=result_bz)

#print(np.array_equal(X, result))
toc = time.time()
print(str(toc - tic))
