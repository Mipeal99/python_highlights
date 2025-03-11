# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:36:22 2019

@author: uo265781
"""

import numpy as np

def simpson (a, b, m, f):
    '''m numero de subintervalos'''
    n = 2*m+1 #numero de nodos
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
f = lambda x: 0.5*(x**2)*np.exp(-x)

def adaptacion(f,a,b,tol):
    '''método adaptativo para calcular integral con simpson'''
    alpha=a #setup del intervalo inicial
    beta=b
    I=0
    while alpha<b:
         #estimacion de la integral
        I1=simpson(alpha,beta,1,f) #con1 intervalo
        I2=simpson(alpha,beta,2,f) #con2 intervalos
        error=1/15.*abs(I1-I2)
        if error<tol*(beta-alpha)/(b-a): #si la integral es buena
            I+=I2 #nos añade el cacho bueno de la integral
            alpha=beta #nuevo punto de inicio donde acabó el bueno
            beta=b #el resto
        else:
            beta=(alpha+beta)/2. #cogemos un intervalo más pequeño
    return I

x=adaptacion(f,a,b,1e-15)
print('adaptacion', x)
print('simpson', simpson(a,b,2000,f))
print(1-x)
print(13/np.e**4)