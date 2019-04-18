import sys
import ROOT
import scipy as sp
import scipy.special
import numpy as np
import h5py
import itertools
import multiprocessing
from math import pi
from itertools import product
import time
np.set_printoptions(threshold=sys.maxsize)
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

# We are simulating magnetic field in a chamber of length 10 m (Z0=-5 to Z0=+5)
# The cross section of the chamber is a square of dim 10m x 10m. Inscribed in itwe have many circular coils of negligible dimension that carries current. The radius of these coils is a=5m. The coils are dz=1 cm apart, so there are 1000 coils spanning 10 m length. The magnetic field inside this chamber will be generated in Cartesian coordinates (x, y, z).

MMT_dim	= int(sys.argv[1])	# m
mu_0	= 4*pi*1.e-7		# Tm/A
a	= []*(100*MMT_dim+1)	# m
I	= []*(100*MMT_dim+1)	# A
C	= []*(100*MMT_dim+1)	# (Tm/A * A)
'''
for zi in range(100*MMT_dim+1):
	a.append(0.5)
	I.append(0.0)
	#print(mu_0*(I[zi]/pi))
	C.append(mu_0*(I[zi]/pi))

I[0]		= 1.e5
I[100*MMT_dim]	= 1.e5
C[0]		= mu_0*(I[0]/pi)
C[100*MMT_dim]  = mu_0*(I[100*MMT_dim]/pi)
'''
#print(str(I[0])+'   '+str(I[100*MMT_dim])+'     '+str(C[0])+'     '+str(C[100*MMT_dim]))


'''
for zi in range(0,51):
	I.append(1.e5-zi*2000.0)
	a.append(0.5)	#1.e-2*(1.0+zi*1.0)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(51,101):
	I.append(0.0+(zi-50)*2000.0)
	a.append(0.5)	#1.e-2*(51.0-(zi-50)*1.0)
	#print(str(zi)+'    '+str(I[zi])+'     '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
'''
'''
for zi in range(0,26):
	a.append(0.25)
	I.append(1.e5)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(26,76):
	a.append(0.50)
	I.append(1.e3)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(76,101):
	a.append(0.25)
	I.append(1.e5)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
'''

for zi in range(0,10*MMT_dim+1):
	a.append(0.10)
	I.append(5.e5)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(10*MMT_dim+1,25*MMT_dim+1):
	a.append(0.75)
	I.append(5.e4)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(25*MMT_dim+1,75*MMT_dim+1):
	a.append(1.20)
	I.append(1.e4)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(75*MMT_dim+1,90*MMT_dim+1):
	a.append(0.75)
	I.append(5.e4)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(90*MMT_dim+1,100*MMT_dim+1):
	a.append(0.10)
	I.append(5.e5)
	#print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))


Bx = []
By = []
Bz = []
X0 =0.0
Y0 =0.0
Z0 =0.0

ArrayBx = np.zeros((100*MMT_dim+1,100*MMT_dim+1,100*MMT_dim+1))
ArrayBy = np.zeros((100*MMT_dim+1,100*MMT_dim+1,100*MMT_dim+1))
ArrayBz = np.zeros((100*MMT_dim+1,100*MMT_dim+1,100*MMT_dim+1))

print('------------xxxxxxxxxxxxx-------------')

XX = range(-50*MMT_dim,50*MMT_dim+1)
YY = range(-50*MMT_dim,50*MMT_dim+1)
ZZ = range(100*MMT_dim+1)

#XX = range(-5,6)
#YY = range(-5,6)
#ZZ = range(0,11)

tic = time.time()

paramlist = list(itertools.product(XX,YY,ZZ))

def calculateField(paramlist):
	xi= paramlist[0]
	yj= paramlist[1]
	zk= paramlist[2]
	L = []
	#print('x '+str(x)+', y '+str(y)+', z '+str(z))
	
	x = 1.e-2*xi
	y = 1.e-2*yj
	
	# loop over n loops
	_bx = 0.0
	_by = 0.0
	_bz = 0.0
	
	for z_n in ZZ:		#iterate over the current loops
		z = 1.e-2*(zk - Z0 - z_n)
		
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
	
	L = [xi, yj, zk, _bx, _by, _bz]
	#L = [yj, zk, _bx, _bz]
	return L

#Generate processes equal to the number of cores
pool = multiprocessing.Pool()

#Distribute the parameter sets evenly across the cores
Ls = []
Ls = pool.map(calculateField,paramlist)

with h5py.File('MagFieldData_list_3D.hdf5','w') as hf:
        hf.create_dataset('B_map',data=Ls)
toc = time.time()
print(str(toc - tic))

