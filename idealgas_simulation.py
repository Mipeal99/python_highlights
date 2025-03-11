# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

boxsize=10
v=1 #modulo velocidad
n1=40
n2=10
n=n1+n2 #numero particulas
m1=1 #masa particulas
m2=50
m=np.concatenate((m1*np.ones(n1),m2*np.ones(n2))) #array de masas
deltat=0.01 #intervalo temporal
npasos=2000 #numero pasos
npause=10 #numero de pasos hasta la pausa
tpause=0.01 # duracion de la pausa
rp=0.15 #radio particulas
kb=0.01 #constante de Boltzmann
T1=20 #temperatura inicial
T2=200
ppause=10 #numero de pasos para calcular presiones



r0=np.random.rand(2,n)*boxsize #posición inicial
theta0=2*np.pi*np.random.rand(n)
#ajustada para que no spawnee nada dentro de la caja y se vuelvan locos los condicionales

e0=np.concatenate((np.random.exponential(kb*T1,(2,n1)),np.random.exponential(kb*T2,(2,n2))),axis=1) #energias iniciales
v0=np.sqrt(2*e0/m) #velocidad iniciañ
#velocidad=np.array([v0*np.cos(theta0),v0*np.sin(theta0)]) #velocidades iniciales, normalizadas a v**2
velocidad=np.array([v0[0,:]*np.cos(theta0),v0[1,:]*np.sin(theta0)])
tstart=timer()

fig=plt.figure()
pligera,=plt.plot(r0[0,:n1],r0[1,:n1],'bo')
ppesada,=plt.plot(r0[0,n1:],r0[1,n1:],'ro')
plt.ylim(0,boxsize) #limites para la caja
plt.xlim(0,boxsize)
t=np.zeros(npasos) #array de tiempos
P=np.zeros(npasos)#array de presiones
pbueno=np.array([0]) #array para presiones teniendo en cuenta ppause
tbueno=np.array([0]) #tiempos con ppause
for i in range(npasos):
    
    r0+=velocidad*deltat #cambio de posicion
    t[i]=i*deltat
#    if i%ppause==0:
#        
#        tbueno=np.append(tbueno,i*deltat) #le meto el tiempo que sea
    ptiempo=0
    for k in range(n):
        for j in range(k+1,n):
            if np.linalg.norm(r0[:,k]-r0[:,j])<2*rp:  #choque
                
                uc=(r0[:,k]-r0[:,j])/np.linalg.norm(r0[:,k]-r0[:,j])#vector unitario en direccion choque
                up=np.array([-uc[1],uc[0]]) #vector unitario perpendicular
                
                vc1=np.dot(velocidad[:,k],uc) #velocidad direccion choque de la particula k 
                vc2=np.dot(velocidad[:,j],uc) # lo mismo con la j
                if vc2-vc1<0:
                    continue
                vp1=np.dot(velocidad[:,k],up) #velocidad direccion perp choque de la particula k 
                vp2=np.dot(velocidad[:,j],up) # lo mismo con la j
                vc1prima=(vc1*(m[k]-m[j])+2*vc2*m[j])/(m[k]+m[j])
                vc2prima=(vc2*(m[j]-m[k])+2*vc1*m[k])/(m[k]+m[j])
                velocidad[:,k]=vc1prima*uc+vp1*up
                velocidad[:,j]=vc2prima*uc+vp2*up
        #v'1c= (v1c(m1-m2)+2*v2c*m2)/(m1+m2)
#v2'c= (v2c(m1-m2)+2*vc2*m2)/(m1+m2)
        
        if (r0[0,k]>boxsize and velocidad[0,k]>0) or (r0[0,k]<0.1 and velocidad[0,k]<0): #rebotes con la caja
            ptiempo+=2*m[k]*2*abs(velocidad[0,k])/deltat/boxsize**2
            velocidad[0,k]=-velocidad[0,k]
            
             #presion
        if (r0[1,k]>boxsize and velocidad[1,k]>0) or (r0[1,k]<0.1 and velocidad[1,k]<0):
            ptiempo+=2*m[k]*2*abs(velocidad[1,k])/deltat/boxsize**2
            
            velocidad[1,k]=-velocidad[1,k]
            
    
    
    P[i]=ptiempo
#    if i%ppause==0:
#        pbueno=np.append(pbueno,np.mean(P[i-ppause:i]))
    
    
    if i%npause==0: #limitamos el numero de redibujos
        plt.pause(tpause) #pausa para que lo haga progresivamente
        pligera.set_data(r0[0,:n1],r0[1,:n1]) #cambiar coordenadas, mas eficiente que plots
        ppesada.set_data(r0[0,n1:],r0[1,n1:])
    plt.show()
#calculo energía
E=np.zeros(n)
for i in range(n):
    
    E[i]=0.5*m[i]*np.dot(velocidad[:,i],velocidad[:,i])
    
tfin=timer()
tejec=tfin-tstart #tiempo de ejecucion


EM=np.mean(E) #energia media
T=EM/kb #temperatura

#fig2=plt.figure()
#plt.hist(E,bins=10,range=(0,2*EM))
#plt.show()


fig3=plt.figure()
plt.plot(tbueno,pbueno)
plt.title("Presion")

pmedia=np.mean(P)
nkbt=n*kb*T
PA=pmedia*boxsize**2 #presion normal
PVDW=pmedia*(boxsize**2-np.pi*rp**2) #van der waals

print("La energía media es ",np.mean(E))
print("Temperatura: ",T)
print("tiempo en ejecutar : ", tejec,"segundos")
print("PA=",PA, " nkbt=",nkbt,"cociente= ",PA/nkbt)
print("P(A-abolas)=",PVDW," nkbt=",nkbt,"cociente= ",PVDW/nkbt)

#extra: meter 2 tipos de particulas: n1=50, n2=50, T1=0, T2=200, m1=1, m2=50

#representar las temperaturas y ver si caen hacia 100

