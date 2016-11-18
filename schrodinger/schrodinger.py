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
#the defalut function is sin(x)
y=np.sin(2*np.pi*time/period)
#y=0.5*(3*time**2-1)
#resol (resolution) is the basis function set length (the real length is resol+1)
resol=10
#vx is potential function, use constant here
vx=1.0
#c is a constant
c=1


#This function is used to calculate input wave function coefficients based on basis functions  
def wave_cf(x, y, basis, resol, T):
    #calculate coefficient for Fourier series 
    #basis function is chosen as exp(i*2Pi*k*x), where k=0, 1, 2, 3, ...
    if basis==1:
      #initialize the list 
      psi_cf=np.zeros(resol+1, dtype=np.complex64)
      for i in range(0, resol+1):
         cntemp= 2*y*np.exp(-1j*2*i*np.pi*x/T)
         psi_cf[i]=cntemp.sum()/cntemp.size  
      return psi_cf
    #Legendre branch, use built in function to calulated the coefficients
    else:
      return np.polynomial.legendre.legfit(x, y, resol)
      
#This function is used to calculate H operator matrix for Fourier series
#this matrix is used in project phase 1, to evalute a input function by the opeartor
def oper_matrix(resol, T):
     operator=np.zeros((resol+1, resol+1),dtype=np.complex64)
     for i in range(0,resol+1):
         operator[i][i]=c*4*np.pi**2*i**2/period**2+vx
     return operator
     
#This function is used to reconstruct the function 
#from basis set coefficients (only used in Fourier series)
#input is indpendent varaible x, resolotion resol, and coefficient in Fourrier space coeff.
def f(x, resol, coeff):      
        f = np.array([coeff[i]*np.exp(1j*2*i*np.pi*x/period) for i in range(0,resol+1)])
        #return the funtion value at position x
        return f.sum()
        
#This function is used to calculate coefficients/vector evaluated by the operator
#H_i=h_ij * b_j
def after_cf(x, y, basis, resol, T, vx):
    #Fourier branch
    if basis==1:
        temp=0j
        #initialize the coefficient 
        psi_after_cf=np.zeros(resol+1, dtype=np.complex64)
        
        for i in range(1,resol+1):
             for j in range(1,resol+1):  
                 #operator matrix
                 temp_matrix=oper_matrix(resol, T)
                 #operator matrix * input wave function vector
                 temp=temp+temp_matrix[i][j]*wave_cf(x, y, basis, resol, T)[j]
             #output function vector
             psi_after_cf[i]=temp
             temp=0j
    #Legendre branch
    else:
        psi_after_cf=np.zeros(resol+1, dtype=np.float64)
        #input wave function vector
        ctemp=wave_cf(x, y, basis, resol, T)
        #apply Laplacian on the vector
        psi_after_cf[0:resol-1]=-c*np.polynomial.legendre.legder(ctemp, m=2, scl=1, axis=0)
        #apply potential energy on the vector
        psi_after_cf=psi_after_cf+vx*ctemp
    #return the coeffcient vector
    return psi_after_cf
    
    
#calculate Hamilton matrix H <psi_i| H |psi_j>, for finding eigenvalue and eigenvector,
#where H * C = lambda * S *C
#use varitional method on http://www.physics.metu.edu.tr/~hande/teaching/741-lectures/lecture-01.pdf
def hamilton_matrix(basis, resol, T):
    #fourier serises 
    if basis==1:
      #initilize the matrix
      hamilton=np.zeros((resol+1, resol+1),dtype=np.complex64)
      #this is a diagnoal matrix, just need to calculated h_ii
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
        #loop a
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
    #find the smallest eigenvalue and eigenvector
    #the smallest eigenvalue state
    if state==0:
         for i in range(0, len(w)):
             if w[i]==min(w):
                 #find the index of smallest eigenvalue
                 index=i
   #delete this part
    #if state==1:
     #find the second smallest eigenvalue
     #    temp=1000000
     #    for i in range(0, len(w)):
     #        if  w[i]<temp:
     #            if w[i]!=min(w):                     
     #               index=i
     #               temp=w[i]
    return w[index], v[:,index]

#main loop for phase 1, show the reconstruct function using two basis functions 

#this loop show the original function and the function after operator 
#basis type 1 for Fourier, 2 for Legendre
for basis in range(1,3):
    #Fourier branch
    if basis==1:   
        #calcuate the coeffcient vector
        psi_after_cf=after_cf(time, y, basis, resol,period, vx)
        #restore the function from basis functions
        y2 = np.array([f(t,resol, psi_after_cf).real for t in time])
        title='Fourier'
    #Legendre branch
    else:
         #calcuate the coeffcient vector
        psi_after_cf=after_cf(time, y, basis, resol,period,vx)
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
    #output the smallest energy state state==0
    w, v=eigen(norm_matrix(basis, resol, period), hamilton_matrix(basis, resol, period), resol, 0)
    #Fourier branch
    if basis==1:     
        #reconstruct the function from basis functions
        y2 = np.array([f(t,resol, v).real for t in time])
        title='Fourier'
    #Legendre branch
    else:
        #reconstruct the function from basis functions
        y2=np.zeros(resol+1, dtype=np.float64)
        y2=np.polynomial.legendre.legval(time, v.real, tensor=True)
        title='Legendre'
    #plot ground state
    axes = plt.gca()  
    axes.set_ylim([-1.2,1.2])
    plt.plot(time, y2, label='ground state')
    pylab.legend(loc='best')
    plt.title(title)
    plt.show()


