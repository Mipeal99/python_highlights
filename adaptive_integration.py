# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:36:22 2019

@author: uo265781
"""

import numpy as np

def simpson (a, b, m, f):
    '''Simpson's rule for numerical integration of the function f on the interval (a,b), m is the number of subintervals'''
    n = 2*m+1 #number of nodes
    x = np.linspace(a, b, n)
    h = (b-a)/(m)
    inte = f(x[0]) + f(x[-1])
    for i in range (1, m):
        inte = inte + 2*f(x[2*i])
    for i in range (1,m+1):
        inte = inte + 4*f(x[2*i-1])
    return h/6*inte

a = 0
b = 4
f = lambda x: 0.5*(x**2)*np.exp(-x) #test function

def adaptacion(f,a,b,tol):
    '''adaptative step method to calculate an integral with Simpson's rule'''
    alpha=a #initial interval setup
    beta=b
    I=0 #the integral's value will be stored here and summed over the subintervals
    while alpha<b:
         #estimation of the integral
        I1=simpson(alpha,beta,1,f) #Simpson with 1 subinterval
        I2=simpson(alpha,beta,2,f) #Simpson with 1 subintervals
        error=1/15.*abs(I1-I2) #calculates the difference between both methods
        if error<tol*(beta-alpha)/(b-a): #if the error is small enough i.e. the approximation is good
            I+=I2 #adds the good chunk to the integral
            alpha=beta #sets up a new interval
            beta=b 
        else: #if the approximation is bad, it picks a smaller interval
            beta=(alpha+beta)/2. 
    return I

x=adaptacion(f,a,b,1e-15) #test run
print('adaptacion', x)
print('simpson', simpson(a,b,2000,f))
print(1-x)

print(13/np.e**4)
