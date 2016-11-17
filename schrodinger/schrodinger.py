# -*- coding: utf-8 -*-
#this is a schrodinger equation project
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import pylab
import os.path
parent = os.path.abspath(os.path.join(os.path.dirname(__file__),'.'))
sys.path.append(parent)

#period. now assume domain is (-period/2, period/2) e,g.(-1,1) in this case
period=2
#generate input wave function
time=np.linspace(-period/2,period/2,100)
y=np.sin(2*np.pi*time/period)
#y=0.5*(3*time**2-1)
#resol is the basis function set length (the real length is resol+1)
resol=20
#vx is potential function, use constant here
vx=2.0
#c is a constant
c=1

#calculate input wave function coefficients    
#store cn(i) in a list
def wave_cf(x, y, basis, resol, T):
    # calculate coefficient for Fourier series 
    #basis function is chosen as exp(i*2Pi*k*x), where k=0, 1, 2, 3, ...
    if basis==1:
      psi_cf=np.zeros(resol+1, dtype=np.complex64)
      for i in range(1, resol+1):
         cntemp= y*np.exp(-1j*2*i*np.pi*x/T)
         psi_cf[i]=cntemp.sum()/cntemp.size  
      # the first one is coefficient for constant basis function 
      psi_cf[0]=y.sum()/y.size
      return psi_cf
    #Legendre branch, use built in function
    else:
      return np.polynomial.legendre.legfit(x, y, resol)
      
#calculate H operator matrix elements h_ij for Fourier series
def oper_matrix(resol, T):
     operator=np.zeros((resol+1, resol+1),dtype=np.complex64)
     for i in range(0,resol+1):
         operator[i][i]=c*4*np.pi**2*i**2/period**2+vx
     return operator
     
#reconstruct the function from basis set coefficients (only used in Fourier series)
def f(x, resol, coeff):      
        f = np.array([2*coeff[i]*np.exp(1j*2*i*np.pi*x/period) for i in range(0,resol+1)])
        return f.sum()
        
#calculate coefficient after operator
def after_cf(x, y, basis, resol, T):
    #Fourier branch
    if basis==1:
        temp=0j
        #initial the coefficient after operator
        psi_after_cf=np.zeros(resol+1, dtype=np.complex64)
        #calculate H matrix elements h_ij
        for i in range(1,resol+1):
             for j in range(1,resol+1):                  
                 temp_matrix=oper_matrix(resol, T)
                 temp=temp+temp_matrix[i][j]*wave_cf(x, y, basis, resol, T)[j]
             psi_after_cf[i]=temp
             temp=0j
    #Legendre branch
    else:
        psi_after_cf=np.zeros(resol+1, dtype=np.float64)
        ctemp=wave_cf(x, y, basis, resol, T)
        #calculate Laplacian on coefficients
        psi_after_cf[0:resol-1]=-c*np.polynomial.legendre.legder(ctemp, m=2, scl=1, axis=0)
        #calculate potential energy on coefficients
        psi_after_cf=psi_after_cf+vx*ctemp
    return psi_after_cf
    
    
#calculate Hamilton matrix H <psi_i| H |psi_j>, for finding eigenvalue and eigenvector,
#where H * C = lambda * S *C
#use varitional method on http://www.physics.metu.edu.tr/~hande/teaching/741-lectures/lecture-01.pdf
def hamilton_matrix(basis, resol, T):
    #fourier serises 
    if basis==1:
      #initilize the matrix
      hamilton=np.zeros((resol+1, resol+1),dtype=np.complex64)
      #this is a diagnoal matrix
      for i in range(0,resol+1):
        hamilton[i][i]=T*(c*4*np.pi**2*i**2/T**2+vx)
    #calcualte hamilton matrix for lengredre polynominal
    else:
      #initilize the matrix
      hamilton=np.zeros((resol+1, resol+1),dtype=np.float64)
      #loop over i
      for i in range(0, resol+1):
        #ctemp1 corresponds to psi_i
        ctempi=np.zeros(resol+1, dtype=np.float64)
        ctempi[i]=1
        #loop over j
        #ctemp2 coressponds to psi_j
        for j in range(0, resol+1):
            ctempj=np.zeros(resol+1, dtype=np.float64)
            ctempj[j]=1
            #laplacian on psi_j,
            ctempj[0:resol-1]=-c*np.polynomial.legendre.legder(ctempj, m=2, scl=1, axis=0) 
            #the list length will be shorter after derivertive, append 0 to the  coeffcient list
            ctempj[resol-1]=0
            ctempj[resol]=0
            #multiply laplacian psi_j with psi_i
            ctemp_sq=np.polynomial.legendre.legmul(ctempi, ctempj)
            #calculate coefficients plus coefficients of potential energy vx times psi_i times psi_j
            ctemp_sq=ctemp_sq+vx*np.polynomial.legendre.legmul(ctempi, ctempj)
            #calculate integral
            ctemp_int=np.polynomial.legendre.legint(ctemp_sq)
            #calculate the integral based on values at each point
            hamilton[i][j]=np.polynomial.legendre.legval(T/2, ctemp_int, tensor=True)-np.polynomial.legendre.legval(-T/2, ctemp_int, tensor=True)
    return hamilton
        
#define S matrix (call it norm matrix here) in varitional method, where H * C = lambda * S *C
#use varitoal method of http://www.physics.metu.edu.tr/~hande/teaching/741-lectures/lecture-01.pdf
def norm_matrix(basis, resol, T):
    #matrix for Fourier
    if basis==1:
        norm=np.zeros((resol+1, resol+1),dtype=np.complex64)
        for i in range(0,resol+1):
            norm[i][i]=T
    #for lengendre
    else:
        norm=np.zeros((resol+1, resol+1),dtype=np.float64)
        for i in range(0, resol+1):
            #set the coefficient for nth order polynominal, eg, ctemp=[0, 1, 0, 0, ....] for 2nd order
            ctempi=np.zeros(resol+1, dtype=np.float64)
            ctempi[i]=1
            #multiply nth order legendre with itself psi_i times psi_i
            ctemp_sq=np.polynomial.legendre.legmul(ctempi, ctempi)
            #evluate the function values at each grid point
            norm_temp=np.polynomial.legendre.legval(time, ctemp_sq, tensor=True)
            #calculate the integral
            norm[i][i]=T*norm_temp.sum()/norm_temp.size
    return norm
    
#calculate eigenvalue and eigenvector of S^-1 * H
def eigen(s, h, resol, state):
    w, v=np.linalg.eig(np.dot(np.linalg.inv(s), h))
    #find the ground state
    #find the two smallest eigenvalue and eigenvector
    #the smallest eigenvalue state
    if state==0:
         for i in range(0, len(w)):
             if w[i]==min(w):
                 index=i
    if state==1:
     #find the second smallest eigenvalue
         temp=1000000
         for i in range(0, len(w)):
             if  w[i]<temp:
                 if w[i]!=min(w):                     
                    index=i
                    temp=w[i]
    return w[index], v[:,index]

#main loop for phase 1, show the reconstruct function using two basis functions 

#this loop show the original function and the function after operator 
#basis type 1 for Fourier, 2 for Legendre
for basis in range(1,3):
    #Fourier branch
    if basis==1:   
        #calcuate the coeffcient vector
        psi_after_cf=after_cf(time, y, basis, resol,period)
        #restore the function from basis functions
        y2 = np.array([f(t,resol, psi_after_cf).real for t in time])
        title='Fourier'
    #Legendre branch
    else:
         #calcuate the coeffcient vector
        psi_after_cf=after_cf(time, y, basis, resol,period)
        y2=np.zeros(resol+1, dtype=np.float64)
        #restore the function from basis functions
        y2 = np.real(np.polynomial.legendre.legval(time, psi_after_cf))
        title='Legendre'
    #output original function and  function  operator  
    plt.plot(time, y, label='original function')
    plt.plot(time, y2, label='operator on the function')
    pylab.legend(loc='best')
    plt.title(title)
    plt.show()

#main loop for phase 2, this loop is used to find the ground state
#basis type 1 for Fourier, 2 for Legendre
for basis in range(1,3):
    #find the eigenvalue and eigenvector
    #output two smallest energy state state==0, state==1
  for state in range(2):
    w, v=eigen(norm_matrix(basis, resol, period), hamilton_matrix(basis, resol, period), resol, state)
    #Fourier branch
    if basis==1:     
        #reconstruct the function from basis functions
        y2 = np.array([f(t,resol, v).real for t in time])
        title='Fourier'
    #Legendre branch
    else:
        #reconstruct the function from basis functions
        y2=np.zeros(resol+1, dtype=np.float64)
        y2=np.polynomial.legendre.legval(time, v, tensor=True)
        title='Legendre'
    #plot ground state
    if state==0:
        plt.plot(time, y2, label='ground state')
    if state==1:
        plt.plot(time, y2, label='second smallest state')
    pylab.legend(loc='best')
    plt.title(title)
    plt.show()


