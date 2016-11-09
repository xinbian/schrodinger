# -*- coding: utf-8 -*-
#this is a schrodinger euqation project


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os.path
parent = os.path.abspath(os.path.join(os.path.dirname(__file__),'.'))
sys.path.append(parent)

#first function, calculate coefficients for fourier series
#the basis set function {constant, cos kx, sin kx}, where k=1, 2, 3, ...

#period. now assume domain is (-period/2, period/2) e,g.(-1,1) in this case
period=2
#generate wave function
time=np.linspace(-period/2,period/2,100)
y=np.sin(2*np.pi*time/period)
#resol is the basis function set length (need to plus 1 for the real length)
resol=50
#vx is potential function, use constant here
vx=2
#c is a constant
c=1
#basis function type 1 for fourrier, 2 for lengendre
basis=2

#store coeffcient and basis function
#initial wave function coefficients
psi_cf=np.zeros(resol+1, dtype=np.complex64)
#wave function after opearator cofficients
#psi_after_cf=np.zeros(resol+1, dtype=np.complex64)
#initial opearator matrix
hamilton=np.zeros((resol+1, resol+1),dtype=np.complex64)
 
    
#store cn(i) in a list, the first one is for constant basis function 
def wave_cf(x, y, basis, resol):
    # calculate coefficient for fourrier serises 
    if basis==1:
      psi_cf=np.zeros(resol+1, dtype=np.complex64)
      for i in range(1, resol+1):
         cntemp= y*np.exp(-1j*2*i*np.pi*time/period)
         psi_cf[i]=cntemp.sum()/cntemp.size  
      # the first one is coefficient for constant basis function 
      psi_cf[0]=y.sum()/y.size
      return psi_cf
    #legendre branch, use built in function
    else:
      return np.polynomial.legendre.legfit(x, y, resol)

#calculate hamilton matrix
def hamilton_matrix(basis, resol):
  for i in range(1,resol+1):
        hamilton[i][i]=c*4*np.pi**2*i**2/period**2+vx
        hamilton[0][0]=vx
  return hamilton

#calculate coeffcient after opearator
def after_cf(basis, resol, x, y):
    #fourrier branch
    if basis==1:
        temp=0j
        #initial the coeffcient after opearator
        psi_after_cf=np.zeros(resol+1, dtype=np.complex64)
        for i in range(1,resol+1):
             for j in range(1,resol+1):  
                 temp=temp+hamilton_matrix(basis, resol)[i][j]*wave_cf(x, y, basis, resol)[j]
             psi_after_cf[i]=temp
             temp=0j
        return psi_after_cf
    #legendre branch
    else:
        ctemp=wave_cf(x, y, basis, resol)
        psi_after_cf=-c*np.polynomial.legendre.legder(ctemp, m=2, scl=1, axis=0)
        psi_after_cf=psi_after_cf+vx*ctemp[0:resol-1]
        return psi_after_cf
            

#reconsttuct the function from basis set coeffcients (only used in fourrier series)
def f(x, Nh, coeff):      
        f = np.array([2*coeff[i]*np.exp(1j*2*i*np.pi*x/period) for i in range(1,Nh+1)])
        return f.sum()
if basis==1:    
    ctemp=wave_cf(time, y, basis, resol)
    psi_after_cf=after_cf(basis, resol, time, y)
    y2 = np.array([f(t,resol, psi_after_cf).real for t in time])+ctemp[0]
else:
    y2=np.polynomial.legendre.legval(time, after_cf(basis, resol, time, y), tensor=True)

plt.plot(time, y)
plt.plot(time, y2)