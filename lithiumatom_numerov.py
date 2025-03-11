# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci

#Parámetros del problema
Z=3 #Número de protones del núcleo
npuntos=2000 #Número puntos en x
infty=30. #Dónde está el infinito

a=Z**2+0.2 #Valor inicial de la energía
da=0.001*a #salto a la hora de hacer el shooting
tol=1e-9 #tolerancia del shooting
nit=70 #iteraciones del shooting
u=np.linspace(infty,1e-6, npuntos) #Array de puntos
h=abs(u[1]-u[0]) #Paso en u


def V(r,n,l,E):
    '''f(x) para el átomo de litio'''
    
    v=((r<1)*1)*(1/r**2*l*(l+1)-2.*3/r+E+4) #Potencial interior
    v+=((r>1)*1)*(1/r**2*l*(l+1)-2./r+E) #Potencial apantallado
    return v

#plt.figure()
#plt.plot(u,V(u,1,0,a))
#plt.xlabel("u",fontsize=14)
#plt.ylabel("f(u)",fontsize=14)
#plt.title("f(u) para el átomo de Li",fontsize=14)

def lacosa(a,npa,n,l):
    p0=np.array([0,h]) #Condiciones iniciales
    def ecsro(a,npa):
        '''Método de Numerov'''

        p=np.array(p0)
        ph=p0*(1-1/12.*h**2*V(u[0:2],n,l,a))
        ph=np.concatenate((ph,np.zeros(npa-2)))
        p=np.concatenate((p,np.zeros(npa-2)))
        for i in range(2,npa):
            ph[i]=2*ph[i-1]-ph[i-2]+h**2*V(u[i-1],n,l,a)*p[i-1]    
            p[i]=ph[i]/(1-1/12.*h**2*V(u[i],n,l,a))
        return p
    
    p= ecsro(a,npa) #Solución inicial
    a0=a 
    da=abs(0.01*a)
    i=1
    while abs(p[-1])>tol and i<400: #Si no es cero en infinito
        a=a-da
        p1=ecsro(a,npa) 
        if p1[-1]*p[-1]<0: #Si hay raiz
            a0+=(i-1)*da
            print("Autovalor encontrado: ","n=",n," l=",l)
            break #Nos vamos al bucle de bisección
        else: #cambiamos el a
            a=a-da       
        i+=1

    i=0 #Reinicio el contador
    p2=np.array([2*tol])
    while i<nit and abs(p2[-1])>tol: #Bisección
        i+=1
        ce=(a+a0)/2.
        p2=ecsro(ce,npa)
        if p2[-1]*p[-1]<0:
            a=ce    
        else:
            a0=ce      
    return p2, ce

nfun=3 # n máximo
autov=np.array([])
#Plot
for n in range(1,nfun+1):
    plt.figure()
    plt.title("n=%i"%int(n))
    plt.xlabel("u=r/$a_0$",fontsize=14)
    plt.ylabel("$U_{n,l}(u)$",fontsize=14)
    pf1,a=lacosa(a-da,npuntos,n,0)
    plt.plot(u,pf1,label="l=0")
    autov=np.append(autov,a)
    for l in range(1,n):
        pf1,b=lacosa(a+10*da,npuntos,n,l)
        autov=np.append(autov,b)       
        plt.plot(u,pf1,label="l=%i" %int(l))

    plt.legend()
print("Autovalores: ",autov)