# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 01:51:20 2023

@author: Usuario
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp #para resolver las ecuaciones

###Parámetros del problema
B=25 #Campo magnetico de Perseus en microgauss
m=30.98e-18 #Masa del ALP en GeV
g=1.01e-11 #Constante de acoplo en GeV^-1
me=0.511e-3 #Masa del electrón en GeV

EE=np.logspace(-5,-1,1200) #Vector de energías para hacer luego los plots y demás

#Densidad de electrones para la frecuencia de plasma
def ne(r): 
    '''ref 49 del paper de FermiLAT ecuación 4, suelta la densidad en cm^-3 y r tiene que estar en kpc'''
    return 3.9e-2/(1+(r/80)**2)**1.8+4.03e-3/(1+(r/280)**2)**0.87
rtest=np.linspace(0,500,500) #Perseus a efectos prácticos para nosotros tiene un radio de 500kpc, por eso el 500
ne=np.max(ne(rtest)) #valor de la densidad de electrones con la que me quedo
wp=4*np.pi*ne*(1.98e-14)**3/137/me #Frecuencia de plasma (al cuadrado)

# def BB(r):
#     return 25*(ne(r)/ne(0))**0.5
# rr=np.linspace(0,1800,500)
# plt.figure()
# plt.plot(rr,BB(rr))

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
tt=np.linspace(0,z,300) #vector de posiciones
cond_inic=[0+0j,1+0j] #condiciones iniciales (el solve_ivp las pide en formato complejo)
i=0
# probs=np.zeros_like(EE)
# for E in EE:#vamos recorriendo todas las energías
#     #Pongo directamente el caso de Mirizzi que es el que estamos usando para este debug
    
#     def eom(z,y):
#         '''y es un vector de variables y z la variable independiente, en este caso escojo que y[0](z) sea el campo del ALP
#         y que y[1](z) sea el A paralelo del fotón.
#         El solve ivp también está hecho de forma que, si el sistema es y'=f(x,y), lo que le tienes que poner es la f
#         así que hay que despejar las derivadas en las ecuaciones del movimiento
#         Por si acaso, j es la unidad imaginaria'''
#         return [-1j*(E*y[0]+deltaa(E)*y[0]+deltaagamma(E)*y[1]),
#                 -1j*(E*y[1]+deltapl(E)*y[1]+deltaagamma(E)*y[0]-1j*g*(E*y[0]+deltaa(E)*y[0]+deltaagamma(E)*y[1])*y[1])]
    
#     sol=solve_ivp(eom,[0,z],cond_inic,t_eval=tt,method='RK45') #resuelvo las ecuaciones del movimiento
#     phi=sol.y[0,:]
#     p0=np.abs(phi[-1])**2
#     probs[i]=1/3*(1-np.exp(-3*Ndominios*p0/2))
#     i+=1
#saco los campos del ALP y del fotón
# 
# aparalelo=sol.y[1,:]
#Plot del numérico, con los valores absolutos al cuadrado
plt.figure(figsize=(12,9))
#plt.plot(tt,np.abs(phi)**2,label='numeric')
########PLOT
#plt.figure(figsize=(12,9))
# plt.plot(1e3*EE,1-probs,label='numeric')

###Ahora calculo y ploteo la probabilidad de Mirizzi

def P0(E,z):
    '''P0 de Mirizzi, z es el tamaño tipico del dominio en kpc (el inverso de las unidades del delta)'''
    return (deltaagamma(E)*z)**2*(np.sin(deltaosc(E)*z/2)/(deltaosc(E)*z/2))**2

def prob(E,z,Ndominios=500):
    '''probabilidad de la exponencial, apéndice del paper de Mirizzi'''
    return 1/3*(1-np.exp(-3*Ndominios*P0(E,z)/2))


z2=500/150
z3=1
# plt.plot(1e3*EE,1-prob(EE,z,50),label='N=50',alpha=0.5)
# plt.plot(1e3*EE,1-prob(EE,z2,150),label='N=150',alpha=0.7)
# plt.plot(1e3*EE,1-prob(EE,z3,500),label='N=500',alpha=0.9)
plt.plot(1e3*EE,1-prob(EE,z),label='N=500',alpha=0.9)
plt.plot([ecrit,ecrit],[0,2],'r--')
plt.xlabel("E(GeV)",fontsize=26)
plt.ylabel("Photon survival probability",fontsize=26)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
# plt.plot(EE,prob(r,EE,z))
plt.ylim(0.6,1.02)
plt.xscale('log')
plt.legend(fontsize=24)

EE2=np.logspace(-4,0,1200)
plt.figure(figsize=(12,9))
plt.xlabel("E(GeV)",fontsize=26)
plt.ylabel(u"Mixing angle (rad)",fontsize=26)
# plt.ylabel(u'$sin^2(\Delta_{osc}z/2)/(\Delta_{osc}z/2)^2$',fontsize=26)
# plt.ylabel('Oscillation length (kpc)',fontsize=28)
plt.xticks(fontsize=23)
plt.yticks(fontsize=23)
plt.plot(1e3*EE2,0.5*np.arctan(2*deltaagamma(EE2)/(deltapl(EE2)-deltaa(EE2))))
# plt.plot([0.1,1000],[np.pi/4,np.pi/4],'g--')
# plt.plot(1e3*EE2,(np.sin(deltaosc(EE2)*z/2)/(deltaosc(EE2)*z/2))**2)
# plt.plot(1e3*EE2,2*np.pi/deltaosc(EE2))
plt.plot(1e3*EE2,np.pi/4*np.ones_like(EE2))
plt.plot([ecrit,ecrit],[0,8.2],'r--')
# plt.ylim(0,8.3)
plt.ylim(0,0.8)
plt.xscale('log')



