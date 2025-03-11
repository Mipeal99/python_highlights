# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#Parámetros

c = 3.0e+8 #Velocidad de la Luz
wl = 632.8e-9 #Longitud de onda
tp = wl/c #Periodo
Nx = 61 #Puntos en x
Nz = 61 #Puntos en z
kd=30 #Límite entre medios/Posición de las rendijas
Nt = 400 #Pasos de tiempo
l0 = 45 #Posición fuente en x
k0 = 15 #Posición fuente en z
dx = wl/7 #Espaciado de los puntos en x
dz = wl/7 #Espaciado de los puntos en z
dt = (dx+dz)/(4*c) #Espaciado temporal
t = 0. #Tiempo inicial
x = dx*np.linspace(1, Nx, Nx)
z = dz*np.linspace(1,Nz,Nz)
zi = dz/1e-6 #Valor inicial de z (micras)
zf = Nz*dz/1e-6 #Valor final de z
xi = dx/1e-6 #Valor inicial de x
xf = Nx*dx/1e-6 #Valor final de x
x, z = np.meshgrid(x, z) #Mallado (x,z)
t0 = 0.75*tp #Tiempo de retardo
sigma = 0.25*tp #Ancho de la fuente gaussiana
ta = 2.0 #Aceleración del tiempo
erel1=4 #Permitividad medio arriba
erel2=1#Permitividad medio abajo
#para incidencia a 45º, reflexión total si erel1=1, erel2=2
eps0=8.85e-12 #Permitividad del vacio
cond1=0 #Conductividades
cond2=0
ctc1 = dt*cond1/(2*erel1*eps0) #Coeficientes teniendo en cuenta
ctc2=dt*cond2/(2*erel1*eps0)  #la conductividad en ambos medios

#Definición de los campos

ey = np.zeros((Nx, Nz)) #Campo eléctrico (y)
hx = np.zeros((Nx,Nz)) #Campo magnético (x)
hz = np.zeros((Nx,Nz)) #Campo magnético (z)

#Preparo la representación
fig = plt.figure(figsize=(11,9))
ax = fig.add_subplot(1,1,1, projection= '3d') #Creamos figura 3d
ax.set_xlim([xi, xf])
ax.set_ylim([zi, zf])
ax.set_zlim([-1.1, 1.1])
ax.set_xlabel('x (micras)')
ax.set_ylabel('z (micras)')
ax.set_zlabel('Ey(z,t)')
ax.view_init(90, 270) #Punto de vista del plot
repre=ax.plot_surface(x/1e-6, z/1e-6, ey.T, cmap = cm.coolwarm, antialiased=False)

#Arrays para las condiciones de contorno
ey0x = np.zeros(Nx)
ey0z = np.zeros(Nz)
eyNx = np.zeros(Nx)
eyNz = np.zeros(Nz)
#Bucles
for n in range (1, Nt+1):
    
    #Condiciones de contorno
	ey[:,0] = ey0x
	ey0x = ey[:,1]
	ey[:,-1] = eyNx
	eyNx = ey[:,-2]
	ey[0,:] = ey0z
	ey0z = ey[1,:]
	ey[-1,:] = eyNz
	eyNz = ey[-2,:]
	
	#Condiciones de contorno en las esquinas
	ey[0,0] = ey[1,1]
	ey[0,-1] = ey[1,-2]
	ey[-1,0] = ey[-2,1]
	ey[-1,-1] = ey[-2,-2]
    
	#Modifico el campo eléctrico
	for l in range (1, Nx-1):
		for k in range (1, Nz-1):
			if k>kd: #Si estamos en el medio de arriba
				ca=(1-ctc1)/(1+ctc1)
				cb=0.5/erel1*(1+ctc1)
				ey[l, k] = ca*ey[l, k] + cb*(hx[l,k]-hx[l,k-1]-hz[l,k] + hz[l-1,k])

			else: #Si estamos en el medio de abajo
				ca=(1-ctc2)/(1+ctc2)
				cb=0.5/erel2*(1+ctc2)
				ey[l, k] = ca*ey[l, k] + cb*(hx[l,k]-hx[l,k-1]-hz[l,k] + hz[l-1,k])
	#Elección de la fuente		
	condi = 'blando'
	tipo = 'plana'
	if condi =='blando':
		condin = 1.
	elif condi == 'duro':
		condin = 0.
    #Introducción de la fuente
	if tipo == 'no plana':
		ey[l0,k0] = condin*ey[l0,k0] + np.exp(-0.5*((t-t0)/sigma)**2.0)
	elif tipo == 'plana':
		l1 = 10 #Longitud de la fuente
		for a in range(-l1,l1):
			ey[l0+a,k0+a] = condin*ey[l0+a,k0+a] + np.exp(-0.5*((t-t0)/sigma)**2.0)  #Fuente gaussiana a 45º
#			ey[l0+a,k0+a]=condin*ey[l0+a,k0+a]+0.5*np.sin(0.1*t/dt) #Fuente sinusoidal
        #Fuente sin inclinación
#		ey[l0-l1:l0+l1,k0] = condin*ey[l0-l1:l0+l1,k0] + np.exp(-0.5*((t-t0)/sigma)**2.0)
#		ey[l0-l1:l0+l1,k0] = condin*ey[l0-l1:l0+l1,k0] + np.sin(0.5*(t-t0))
	#Modifico el campo magnético
	for l in range (0, Nx-1):
		for k in range (0, Nz-1):
			hx[l,k] = hx[l,k]+0.5*(ey[l,k+1]-ey[l,k])
			hz[l,k] = hz[l,k]+0.5*(ey[l,k]-ey[l+1,k])
    #Actualizo el plot
	if np.ceil(n/ta)==n/ta:
		repre.remove()
		ax.set_title('n = %i pasos temporales' %n, x=0.11, y=0.99)
		repre=ax.plot_surface(x/1e-6, z/1e-6, ey.T, cmap=cm.coolwarm, antialiased=False)
		plt.pause(0.001)
	t = n*dt
plt.show()
