# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 01:51:20 2023

@author: Usuario
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp #para resolver las ecuaciones
import os

###Parámetros del problema
B=25 #Campo magnetico de Perseus en microgauss
m=3.98e-18 #Masa del ALP en GeV
g=1e-22 #Constante de acoplo en GeV^-1
me=0.511e-3 #Masa del electrón en GeV

EE=np.logspace(-5,1,1800) #Vector de energías para hacer luego los plots y demás


#Densidad de electrones para la frecuencia de plasma
def ne(r): 
    '''ref 49 del paper de FermiLAT ecuación 4, suelta la densidad en cm^-3 y r tiene que estar en kpc'''
    return 3.9e-2/(1+(r/80)**2)**1.8+4.03e-3/(1+(r/280)**2)**0.87
rtest=np.linspace(0,500,500) #Perseus a efectos prácticos para nosotros tiene un radio de 500kpc, por eso el 500
ne=np.max(ne(rtest)) #valor de la densidad de electrones con la que me quedo
wp=4*np.pi*ne*(1.98e-14)**3/137/me #Frecuencia de plasma (al cuadrado)

###DEFINO LOS DELTAS SEGÚN EL PAPER 
###"Hardening of TeV gamma spectrum of AGNs in galaxy clusters by conversions of photons into axion-like particles"
def deltapl(E):
    '''suelta kpc^-1, E va en TeV y ne va en cm^-3'''
    return -1.1e-10/E*ne/1e-3
def deltaagamma(E):
    '''suelta kpc^-1, g va en GeV^-1 y B va en microgauss'''
    return 7.6e-2*g/5e-11*B
def deltaa(E):
    '''suelta kpc^-1, m va en neV y E va en TeV'''
    return -7.8e-5/E*(m*1e18)**2

ecrit=50/100*np.abs((m*1e18)**2-(wp*1e18)**2)/B*5e-11/g #Energia critica en GeV, m y wp tienen que ir en neV, B en microgauss
def deltaosc(E):
    '''delta oscilacion segun la ref 41 de FermiLAT, en terminos de la energia critica, suelta kpc^-1
    '''
    return 2*deltaagamma(E)*np.sqrt(1+(ecrit*1e-3/E)**2)


###RESOLUCIÓN NUMÉRICA

###Parámetros para el tema de la exponencial
Ndominios=300
s=500
z=s/Ndominios #tamaño típico de cada dominio
#como la probabilidad se calcula al final de cada dominio, la integración numérica tiene que ir de 0 a z
# tt=np.linspace(0,z,300) #vector de posiciones
cond_inic=[0+0j,1/np.sqrt(2)+0j,1/np.sqrt(2)+0j] #condiciones iniciales (el solve_ivp las pide en formato complejo)
i=0
probs=np.zeros_like(EE)
pceros=np.zeros_like(EE)
for E in EE:#vamos recorriendo todas las energías
    #DEFINO LAS ECUACIONES DEL MOVIMIENTO
    
    def eom(z,y):
        '''y es un vector de variables y z la variable independiente, en este caso escojo que y[0](z) sea el campo del ALP
        y que y[1](z) sea el A paralelo del fotón.
        El solve ivp también está hecho de forma que, si el sistema es y'=f(x,y), lo que le tienes que poner es la f
        así que hay que despejar las derivadas en las ecuaciones del movimiento
        Por si acaso, j es la unidad imaginaria'''
        partialz=-1j*(E*y[0]+deltaa(E)*y[0]-1e3*g*E/4*(y[1]*y[1]+y[2]*y[2])+9.27e15*g*B**2/4/E)
        return [partialz, -1j*y[1]+1j/(1+g*y[0])*(1e3*g/2*y[0]*y[1]*E-deltapl(E)*y[1]),
                -1j*y[2]+1j/(1+g*y[0])*(1e3*g/2*y[0]*y[2]*E-deltapl(E)*y[2]-g*B/2/E*partialz)]
    
    sol=solve_ivp(eom,[0,z],cond_inic,method='BDF') #resuelvo las ecuaciones del movimiento
    phi=sol.y[0,:]
    p0=np.abs(phi[-1])**2
    pceros[i]=p0
    probs[i]=1/3*(1-np.exp(-3*Ndominios*p0/2))
    i+=1
#saco los campos del ALP y del fotón
# 
# aparalelo=sol.y[1,:]
#Plot del numérico, con los valores absolutos al cuadrado
# plt.figure()
# plt.plot(tt,np.abs(phi)**2,label='numeric')
########PLOT
plt.figure(figsize=(12,9))
plt.plot(1e3*EE,1-probs,label=u'$F_0$ model')

# ruta=os.path.join(r'C:\Users\Usuario\Desktop\pthones\TFM','F0.txt')
# np.save(ruta,1-probs)
###Ahora calculo y ploteo la probabilidad de Mirizzi
###IMPORTANTE: REDEFINO EL g PARA QUE AQUI SALGA EL PLOT AL QUE ESTAMOS ACOSTUMBRADOS

g=1.01e-11

def deltapl(E):
    '''suelta kpc^-1, E va en TeV y ne va en cm^-3'''
    return -1.1e-10/E*ne/1e-3
def deltaagamma(E):
    '''suelta kpc^-1, g va en GeV^-1 y B va en microgauss'''
    return 7.6e-2*g/5e-11*B
def deltaa(E):
    '''suelta kpc^-1, m va en neV y E va en TeV'''
    return -7.8e-5/E*(m*1e18)**2

ecrit=50/100*np.abs((m*1e18)**2-(wp*1e18)**2)/B*5e-11/g #Energia critica en GeV, m y wp tienen que ir en neV, B en microgauss
def deltaosc(E):
    '''delta oscilacion segun la ref 41 de FermiLAT, en terminos de la energia critica, suelta kpc^-1
    '''
    return 2*deltaagamma(E)*np.sqrt(1+(ecrit*1e-3/E)**2)

def P0(E,z):
    '''P0 de Mirizzi, z es el tamaño tipico del dominio en kpc (el inverso de las unidades del delta)'''
    return (deltaagamma(E)*z)**2*(np.sin(deltaosc(E)*z/2)/(deltaosc(E)*z/2))**2

def prob(E,z,Ndominios=500):
    '''probabilidad de la exponencial, apéndice del paper de Mirizzi'''
    return 1/3*(1-np.exp(-3*Ndominios*P0(E,z)/2))

#ploteo Mirizzi
plt.plot(1e3*EE,1-prob(EE,z),label=u'$G_0$ model')
plt.xlabel("E(GeV)",fontsize=28)
plt.ylabel("Photon survival probability",fontsize=28)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
# plt.plot(EE,prob(r,EE,z))
plt.ylim(0.6,1.02)
plt.xscale('log')
plt.legend(fontsize=26)




