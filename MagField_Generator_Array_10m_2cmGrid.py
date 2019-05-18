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

#---------------------------------------------------------
for zi in range(-50*MMT_dim,50*MMT_dim+1):
	a.append(0.0)
	I.append(0.0)
	C.append(0.0)
#---------------------------------------------------------

for zi in range(-50*MMT_dim,-40*MMT_dim,10):
	fi = abs((abs(zi)-500)/10)
	a[zi] = float("%0.2f" % (0.05 + 0.1*fi))#0.05 + 0.1*fi
	I[zi] = 8.e5 -4.e4*fi 
	C[zi] = mu_0*(I[zi]/pi)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))

for zi in range(-40*MMT_dim,40*MMT_dim+1,20):
	a[zi] = 7.50
	I[zi] = 4.e5
	C[zi] = mu_0*(I[zi]/pi)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))

for zi in range(41*MMT_dim,50*MMT_dim+1,10):
	fi = 10-(500-abs(zi))/10
	a[zi] = float("%0.2f" % (1.05 - 0.1*fi))#1.05 - 0.1*fi
	I[zi] = 4.e5 +4.e4*fi
	C[zi] = mu_0*(I[zi]/pi)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))



XX = range(-50*MMT_dim,50*MMT_dim+1,2)
YY = range(-50*MMT_dim,50*MMT_dim+1,2)
ZZ = range(-50*MMT_dim,50*MMT_dim+1,2)

tic = time.time()

paramlist = list(itertools.product(XX,YY,ZZ))

print(str(type(paramlist)))


result_bx = np.ctypeslib.as_ctypes(np.empty((50*MMT_dim+1,50*MMT_dim+1,50*MMT_dim+1)))
shared_array_bx = sharedctypes.RawArray(result_bx._type_, result_bx)

result_by = np.ctypeslib.as_ctypes(np.empty((50*MMT_dim+1,50*MMT_dim+1,50*MMT_dim+1)))
shared_array_by = sharedctypes.RawArray(result_by._type_, result_by)

result_bz = np.ctypeslib.as_ctypes(np.empty((50*MMT_dim+1,50*MMT_dim+1,50*MMT_dim+1)))
shared_array_bz = sharedctypes.RawArray(result_bz._type_, result_bz)

def calculateField(paramlist):
	xi= paramlist[0]
	yj= paramlist[1]
	zk= paramlist[2]
	
	tmp_bx = np.ctypeslib.as_array(shared_array_bx)
	tmp_by = np.ctypeslib.as_array(shared_array_by)
	tmp_bz = np.ctypeslib.as_array(shared_array_bz)
	
	x = float("%0.2f" % (1.e-2*xi))	#1.e-2*xi
	y = float("%0.2f" % (1.e-2*yj))	#1.e-2*yj
	
	_bx = 0.0
	_by = 0.0
	_bz = 0.0
	
	for z_n in ZZ:          #iterate over the current loops
		if a[z_n]==0.0:
			continue
		z = float("%0.2f" % (1.e-2*(zk - z_n))) #1.e-2*(zk - z_n)
		
		if x==0.0 and y==0.0:
			x = 1.e-6
			y = 1.e-6
		#if x==0.0:
		#	x = 1.e-3
		#if y==0.0:
		#	y = 1.e-3
		#if z==0.0:
		#	z = 1.e-3
		
		r       = np.sqrt(x*x + y*y + z*z)
		rho     = np.sqrt(x*x + y*y)
		alpha   = np.sqrt(a[z_n]*a[z_n] + r*r - 2*a[z_n]*rho)
		beta    = np.sqrt(a[z_n]*a[z_n] + r*r + 2*a[z_n]*rho)
		k       = np.sqrt(1.0-((alpha*alpha)/(beta*beta)))
		
		if k>0.95:
			k=0.95
		
		if alpha<0.1:
			continue
		
		_bx += ((C[z_n]*x*z)/(2*alpha*alpha*beta*rho*rho))*((a[z_n]*a[z_n]+r*r)*scipy.special.ellipe(k*k) - alpha*alpha*scipy.special.ellipk(k*k))
		_by += ((C[z_n]*y*z)/(2*alpha*alpha*beta*rho*rho))*((a[z_n]*a[z_n]+r*r)*scipy.special.ellipe(k*k) - alpha*alpha*scipy.special.ellipk(k*k))
		_bz += (C[z_n]/(2*alpha*alpha*beta))*((a[z_n]*a[z_n]-r*r)*scipy.special.ellipe(k*k) + alpha*alpha*scipy.special.ellipk(k*k))
		
		_b = np.sqrt(_bx*_bx + _by*_by + _bz*_bz)
		if _b>15:
			print(str(x)+'     '+str(y)+'     '+str(z)+'     '+str(zk)+'     '+str(z_n)+'     '+str(a[z_n])+'     '+str(_bx)+'     '+str(_by)+'     '+str(_bz)+'     '+str(alpha)+'     '+str(beta)+'     '+str(k))
		
	tmp_bx[xi, yj, zk] = _bx
	tmp_by[xi, yj, zk] = _by
	tmp_bz[xi, yj, zk] = _bz
	#print(str(xi)+'     '+str(yj)+'     '+str(zk)+'     '+str(_bx)+'     '+str(_by)+'     '+str(_bz))

p = Pool()
res = p.map(calculateField, paramlist)
result_bx = np.ctypeslib.as_array(shared_array_bx)
result_by = np.ctypeslib.as_array(shared_array_by)
result_bz = np.ctypeslib.as_array(shared_array_bz)

with h5py.File('MagFieldData_array_2cm_grid_gpu.hdf5','w') as hf:
	hf.create_dataset('Bx_map',data=result_bx)
	hf.create_dataset('By_map',data=result_by)
	hf.create_dataset('Bz_map',data=result_bz)

toc = time.time()
print(str(toc - tic))
