# -*- coding: utf-8 -*-
#this is a schrodinger euqation project


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os.path
import pylab
parent = os.path.abspath(os.path.join(os.path.dirname(__file__),'.'))
sys.path.append(parent)

#first function, calculate coefficients for fourier series
#the basis set function {constant, cos kx, sin kx}, where k=1, 2, 3, ...

#period. now assume domain is (0, period)
period=2*np.pi
time=np.linspace(-period/2,period/2,100)
y=np.sin(time)
#resol is the basis function set length (need to plus1 for the real length if dft)
resol=50
#vx is potential function
vx=2
#c is a constant
c=3
#basis function type
basis=2

#store coeffcient and basis function
#wave function coefficients
psi_cf=np.zeros(resol+1, dtype=np.complex64)
#wave function after opearator cofficients
#psi_after_cf=np.zeros(resol+1, dtype=np.complex64)
#opearator matrix
hamilton=np.zeros((resol+1, resol+1),dtype=np.complex64)


def cn(y, n, basis):
    #fourrier serises branch
    if basis==1:
        c = y*np.exp(-1j*2*n*np.pi*time/period)
        return c.sum()/c.size  
    


#store cn(i) in a list, the first one is for constant coeffcients
def wave_cf(x, y, basis, resol):
    # calculate constant coefficient for fourrier serises 
    if basis==1:
      psi_cf=np.zeros(resol+1, dtype=np.complex64)
      for i in range(1, resol+1):
         psi_cf[i]=cn(y, i, basis)
      psi_cf[0]=y.sum()/y.size
      return psi_cf
     #legendre branch
    else:
      return np.polynomial.legendre.legfit(time, y, resol)

#calculate hamilton matrix
def hamilton_matrix(basis, resol):
  for i in range(1,resol+1):
        hamilton[i][i]=c*4*np.pi**2*i**2/period**2+vx
        hamilton[0][0]=vx
  return hamilton

#calculate coeffcient after opearator
temp=0j
if basis==1:
 psi_after_cf=np.zeros(resol+1, dtype=np.complex64)
 for i in range(1,resol+1):
    for j in range(1,resol+1):
    
           temp=temp+hamilton_matrix(basis, resol)[i][j]*wave_cf(time, y, basis, resol)[j]
    psi_after_cf[i]=temp
    temp=0j
    
if basis==1:
    pass
else:
    ctemp=wave_cf(time, y, basis, resol)
    psi_after_cf=-c*np.polynomial.legendre.legder(ctemp, m=2, scl=1, axis=0)
    psi_after_cf=psi_after_cf+vx*ctemp[0:resol-1]
    
#cder=np.polynomial.legendre.legder(c2, m=2, scl=1, axis=0)
#y2=np.polynomial.legendre.legval(time, cder, tensor=True)

def f(x, Nh, coeff):      
        f = np.array([2*coeff[i]*np.exp(1j*2*i*np.pi*x/period) for i in range(1,Nh+1)])
        return f.sum()

if basis==1:    
    ctemp=wave_cf(time, y, basis, resol)
    y2 = np.array([f(t,resol, psi_after_cf).real for t in time])+ctemp[0]
else:
    y2=np.polynomial.legendre.legval(time, psi_after_cf, tensor=True)

plt.plot(time, y)
plt.plot(time, y2)