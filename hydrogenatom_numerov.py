# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
#Parametros del problema
npuntos=5000 #Número puntos en x
infty=50. #Dónde esta el infinito

a=1.2 #Valor inicial de la energía
da=0.001*a #Salto a la hora de hacer el shooting
tol=1e-9 #Tolerancia del shooting
nit=70 #Iteraciones del shooting
u=np.linspace(infty,1e-6, npuntos) #Array de puntos
h=abs(u[1]-u[0]) #Paso en u

def V(r,n,l,E):
    '''f(x) para el átomo de hidrógeno'''
    return (1/r**2*l*(l+1)-2./r+E)

def lacosa(a,npa,n,l):
    
    p0=np.array([0,h]) #Condiciones iniciales
    def ecsro(a,npa):
        '''Método de Numerov'''
    
        p=np.array(p0)
        ph=p0*(1-1/12.*h**2*V(u[0:2],n,l,a)) #Pasamos a las phi
        ph=np.concatenate((ph,np.zeros(npa-2))) #Relleno con ceros
        p=np.concatenate((p,np.zeros(npa-2)))
        for i in range(2,npa): #Aplico la fórmula
            ph[i]=2*ph[i-1]-ph[i-2]+h**2*V(u[i-1],n,l,a)*p[i-1]    
            p[i]=ph[i]/(1-1/12.*h**2*V(u[i],n,l,a))
        return p 
    
    p= ecsro(a,npa) #Solución inicial
    
    a0=a #Guardo el valor de a
    da=abs(0.01*a) #Redefino el salto para que no de problemas
    i=1 #Contador
    while abs(p[-1])>tol and i<400: #Si no es cero en infinito
        a=a-da #Como las energías son negativas y van acercándose a cero, a disminuye
        p1=ecsro(a,npa) #Calcula función de onda para otra energia
        if p1[-1]*p[-1]<0: #Si hay raiz
            a0+=(i-1)*da
            print("Autovalor encontrado: ","n=",n," l=",l)
            break #Nos vamos al bucle de bisección
        else: #Cambiamos el a
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

nfun=3 # n máximo para sacar funciones (saca nfun figuras con 2nfun-1 plots)
autov=np.array([])
#Ploteo las funciones para los valores de l permitidos
for n in range(1,nfun+1):
    plt.figure()
    plt.title("n=%i"%int(n))
    pf1,a=lacosa(a-da,npuntos,n,0)
    plt.plot(u,pf1,label="l=0")
    autov=np.append(autov,a)
    for l in range(1,n):
        pf1,b=lacosa(a+10*da,npuntos,n,l)
        autov=np.append(autov,b)       
        plt.plot(u,pf1,label="l=%i" %int(l))  
    plt.legend()
enes=np.sqrt(1/autov) #Revierto el cálculo de la energía al número cuántico n
print("n= ",enes)