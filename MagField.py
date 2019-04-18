import sys
import ROOT
import scipy as sp
import scipy.special
import numpy as np
import h5py
from math import pi

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
	a.append(0.25)
	I.append(1.e5)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(10*MMT_dim+1,25*MMT_dim+1):
	a.append(0.75)
	I.append(1.e4)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(25*MMT_dim+1,75*MMT_dim+1):
	a.append(2.00)
	I.append(1.e3)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(75*MMT_dim+1,90*MMT_dim+1):
	a.append(0.75)
	I.append(1.e4)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
	C.append(mu_0*(I[zi]/pi))
for zi in range(90*MMT_dim+1,100*MMT_dim+1):
	a.append(0.25)
	I.append(1.e5)
	print(str(zi)+'    '+str(I[zi])+'    '+str(a[zi]))
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

for xi in range(-50*MMT_dim,50*MMT_dim+1):					#-50*MMT_dim,50*MMT_dim
	for yj in range(-50*MMT_dim,50*MMT_dim+1):				#-50*MMT_dim,50*MMT_dim
		for zk in range(100*MMT_dim+1):
			#print(zk)
			# @point (xi,yi,zk)
			x = 1.e-2*(xi - X0)
			y = 1.e-2*(yj - Y0)
			
			# loop over n loops
			_bx = 0.0
			_by = 0.0
			_bz = 0.0
			for z_n in range(100*MMT_dim+1):
				z = 1.e-2*(zk - Z0 - z_n)
				
				r	= np.sqrt(x*x + y*y + z*z)
				rho	= np.sqrt(x*x + y*y)
				alpha	= np.sqrt(a[z_n]*a[z_n] + r*r - 2*a[z_n]*rho)
				beta	= np.sqrt(a[z_n]*a[z_n] + r*r + 2*a[z_n]*rho)
				k       = np.sqrt(1.0-((alpha*alpha)/(beta*beta)))
				gamma   = x*x - y*y
				
				_bx += ((C[z_n]*x*z)/(2*alpha*alpha*beta*rho*rho))*((a[z_n]*a[z_n]+r*r)*scipy.special.ellipe(k*k) - alpha*alpha*scipy.special.ellipk(k*k))
				_by += ((C[z_n]*y*z)/(2*alpha*alpha*beta*rho*rho))*((a[z_n]*a[z_n]+r*r)*scipy.special.ellipe(k*k) - alpha*alpha*scipy.special.ellipk(k*k))
				_bz += (C[z_n]/(2*alpha*alpha*beta))*((a[z_n]*a[z_n]-r*r)*scipy.special.ellipe(k*k) + alpha*alpha*scipy.special.ellipk(k*k))
			
			if np.isnan(_bx)==True:
				_bx = 1.e-9
			if np.isnan(_by)==True:
				_by = 1.e-9
			if np.isnan(_bz)==True:
				_bz = 1.e-9
			
			ArrayBx[xi,yj,zk]=_bx
			ArrayBy[xi,yj,zk]=_by
			ArrayBz[xi,yj,zk]=_bz
	
with h5py.File('MagFieldData_array_10m.hdf5','w') as hf:
	hf.create_dataset('Bx_map',data=ArrayBx)
	hf.create_dataset('By_map',data=ArrayBy)
	hf.create_dataset('Bz_map',data=ArrayBz)

