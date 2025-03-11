# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 02:45:48 2023

@author: Miguel
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci
import scipy.stats as stats
import time


start_time = time.time()


###Parámetros del problema
B=10 #Campo magnetico de Perseus en microgauss

me=0.511e-3 #Masa del electrón en GeV

masas=np.logspace(-18,-15,10)
masasF=masas[12:]
masas2=np.logspace(-18,-13,200)
gs=np.logspace(-13,-8,90)
gtildes=np.logspace(-20.5,-16,20)
gtildes=gtildes[:12]
g1s=np.logspace(-13,-10,10)
g2s=[1e43,1e44,1e45,1e46]

#g=0.5e-11 #Constante de acoplo del modelo G0 en GeV^-1
#g2=0.5e42 #G2 del modelo G0^2 en GeV^-5
#g1=0.6e-11 #G1 del modelo G0^2 en GeV^-1
#gtilde=2e-19 #g del modelo F0 en GeV^-1
###Parámetros de la power law
E0=0.3 #TeV
Ecut=0.560 #TeV
N=1.54e-9
gamma=2.11




#Bins de energia
nbins=177
bineado=np.logspace(-1,1,nbins+1) #inicios y finales de cada bin en energia (TeV)



# EE=np.logspace(-1,1,500) #Vector de energías

#Densidad de electrones para la frecuencia de plasma
def ne(r): 
    '''ref 49 del paper de FermiLAT ecuación 4, suelta la densidad en cm^-3 y r tiene que estar en kpc'''
    return 3.9e-2/(1+(r/80)**2)**1.8+4.03e-3/(1+(r/280)**2)**0.87
rtest=np.linspace(0,500,500) #Perseus a efectos prácticos para nosotros tiene un radio de 500kpc, por eso el 500
ne=np.max(ne(rtest)) #valor de la densidad de electrones con la que me quedo
wp=4*np.pi*ne*(1.98e-14)**3/137/me #Frecuencia de plasma (al cuadrado)

#Dominios (físicos) en los que el campo magnetico es igual/tiene la misma dirección (apéndice de Mirizzi/sección 3.1 del TFM)
Ndominios=600
s=1865 #extensión total que consideramos de Perseus en kpc (radio virial)
z=s/Ndominios #tamaño de cada dominio

###DEFINO LOS DELTAS SEGÚN EL PAPER de la siguiente linea (comprobado por mí ya) 
###"Hardening of TeV gamma spectrum of AGNs in galaxy clusters by conversions of photons into axion-like particles"
def deltapl(E):
    '''suelta kpc^-1, E va en TeV y ne va en cm^-3'''
    return -1.1e-10/E*ne/1e-3
def deltaagamma(E,g):
    '''suelta kpc^-1, g va en GeV^-1 y B va en microgauss'''
    return 7.6e-2*g/5e-11*B
def deltaa(E,m):
    '''suelta kpc^-1, m va en neV y E va en TeV'''
    return -7.8e-5/E*(m*1e18)**2


def deltaosc(E):
    '''delta oscilacion segun la ref 41 de FermiLAT, en terminos de la energia critica, suelta kpc^-1
    '''
    return 2*deltaagamma(E,g)*np.sqrt(1+(ecrit*1e-3/E)**2)
def P0(E,z):
    '''P0 de Mirizzi, z es el tamaño tipico del dominio en kpc (el inverso de las unidades del delta)'''
    return (deltaagamma(E,g)*z)**2*(np.sin(deltaosc(E)*z/2)/(deltaosc(E)*z/2))**2

def prob_G0(E,z=z,Ndominios=Ndominios):
    '''probabilidad de la exponencial del modelo de G0, apéndice del paper de Mirizzi'''
    return 1-1/3*(1-np.exp(-3*Ndominios*P0(E,z)/2))

def prob_F0(E,z=z,Ndominios=Ndominios):
    '''Programa del modelo de F0 adaptado para sacar las constraints, básicamente copio las ecuaciones del movimiento, resuelvo
    y calculo la probabilidad con la exponencial'''
#    tt=np.linspace(0,z,300) #vector de posiciones para el t_eval si da problemas
    cond_inic=[0+0j,1/np.sqrt(2)+0j,1/np.sqrt(2)+0j] #condiciones iniciales (el solve_ivp las pide en formato complejo)
        #DEFINO LAS ECUACIONES DEL MOVIMIENTO
        
    def eom(z,y):
            '''y es un vector de variables y z la variable independiente, en este caso escojo que y[0](z) sea el campo del ALP
            y que y[1](z) sea el A paralelo del fotón.
            El solve ivp también está hecho de forma que, si el sistema es y'=f(x,y), lo que le tienes que poner es la f
            así que hay que despejar las derivadas en las ecuaciones del movimiento
            Por si acaso, j es la unidad imaginaria'''
            partialz=-1j*(E*y[0]+deltaa(E,m)*y[0]+1e3*gtilde*E/4*(y[1]*y[1]+y[2]*y[2])-9.27e15*gtilde*B**2/4/E)
            return [partialz, -1j*y[1]+1j/(1+gtilde*y[0])*(1e3*gtilde/2*y[0]*y[1]*E-deltapl(E)*y[1]),
                -1j*y[2]+1j/(1+gtilde*y[0])*(1e3*gtilde/2*y[0]*y[2]*E-deltapl(E)*y[2]-gtilde*B/2/E*partialz)]
        
    sol=sci.solve_ivp(eom,[0,z],cond_inic,method='RK45') #resuelvo las ecuaciones del movimiento
    phi=sol.y[0,:] #cojo el campo del ALP
    p0=np.abs(phi[-1])**2 #calculo la probabilidad en un dominio
    prob=1/3*(1-np.exp(-3*Ndominios*p0/2)) #la probabilidad total teniendo en cuenta los N dominios
    return 1-prob

def deltanuevo(E,g2):
    '''el término extra de las ecuaciones de G_0^2'''
    return 1.92e-49*g2*E*B**2

def prob_G0cuadrado(E,z=z,Ndominios=Ndominios):
    '''Programa del modelo de G0^2 adaptado para sacar las constraints, básicamente copio las ecuaciones del movimiento, resuelvo
    y calculo la probabilidad con la exponencial'''
    cond_inic=[0+0j,1+0j]
    def eom(z,y):
        '''y es un vector de variables y z la variable independiente, en este caso escojo que y[0](z) sea el campo del ALP
        y que y[1](z) sea el A paralelo del fotón.
        El solve ivp también está hecho de forma que, si el sistema es y'=f(x,y), lo que le tienes que poner es la f
        así que hay que despejar las derivadas en las ecuaciones del movimiento
        Por si acaso, j es la unidad imaginaria'''
        
        return [-1j*(E*y[0]+deltaa(E,m)*y[0]+deltaagamma(E,g1)*y[1]+deltanuevo(E,g2)*y[1]*y[1]),
                -1j*(E*y[1]+deltapl(E)*y[1]+deltaagamma(E,g1)*y[0]+4*deltanuevo(E,g2)*y[1]*y[0])]
    
    sol=sci.solve_ivp(eom,[0,z],cond_inic,method='BDF') #resuelvo las ecuaciones del movimiento
    phi=sol.y[0,:]
    p0=np.abs(phi[-1])**2
    prob=1/3*(1-np.exp(-3*Ndominios*p0/2))
    return 1-prob



##################################################
####EL PROCESO DE CONSTRAINTS EN SI###############
##################################################


def powerlaw(E,N=N,gamma=gamma,E0=E0,Ecut=Ecut):
    '''para unificarlo todo E va en TeV, el resto son valores del paper'''
    return N*(E/E0)**(gamma)*np.exp(-E/Ecut)

def powerlawconalps_G0(E,N=N,gamma=gamma,E0=E0,Ecut=Ecut):
    '''La misma power law con la photon survival probability para integrar luego'''
    return N*(E/E0)**(gamma)*np.exp(-E/Ecut)*prob_G0(E)


def powerlawconalps_F0(E,N=N,gamma=gamma,E0=E0,Ecut=Ecut):
    '''La misma power law con la photon survival probability para integrar luego'''
    return N*(E/E0)**(gamma)*np.exp(-E/Ecut)*prob_F0(E)

def powerlawconalps_G0cuadrado(E,N=N,gamma=gamma,E0=E0,Ecut=Ecut):
    '''La misma power law con la photon survival probability para integrar luego'''
    return N*(E/E0)**(gamma)*np.exp(-E/Ecut)*prob_G0cuadrado(E)

#creo arrays para almacenar las cosas
cuentas_sin=np.zeros(nbins)
cuentas_conG0=np.zeros(nbins)
cuentas_conF0=np.zeros(nbins)
cuentas_conG02=np.zeros(nbins)

for i in range(nbins):
    #calculo las cuentas esperadas en cada bin según el paper de Fermi/sección 5 del Overleaf
    cuentas_sin[i]=sci.quad(powerlaw,bineado[i],bineado[i+1])[0]
errores=cuentas_sin*0.1 #error del 10% a las cuentas sin
    
#CUANTIL DE LA DISTRIBUCION QUE DEBERIA DE SEGUIR EL CHI CUADRADO
treshold=stats.chi2.ppf(0.05,nbins-1)

rutag=r'C:\Users\Usuario\Desktop\pthones\RESULT_G0.txt'
rutaf=r'C:\Users\Usuario\Desktop\pthones\RESULT_F0_errores20.txt'
rutag2=r'C:\Users\Usuario\Desktop\pthones\RESULT_G02_errores20.txt'

archivog=open(rutag,'a')
archivof=open(rutaf,'a')
archivog2=open(rutag2,'a')

for masa in masas2:
    m=masa
    for cte in gs:
        g=cte
        ecrit=50/100*np.abs((m*1e18)**2-(wp*1e18)**2)/B*5e-11/g #Energia critica en GeV, m y wp tienen que ir en neV, B en microgauss
        for i in range(nbins):
            cuentas_conG0[i]=sci.quad(powerlawconalps_G0,bineado[i],bineado[i+1])[0]
            chicuadrado_G0=np.sum((cuentas_sin-cuentas_conG0)**2/errores**2)
        print('SIGUIENTE MODELO')
        print('g=',g)
        print('masa =', m)
        print('chi cuadrado G0 = ',chicuadrado_G0)
        print('treshold = ',treshold)
        if chicuadrado_G0>treshold:
            print('Excluido')
            archivog.write('g='+str(g)+'\t m='+str(m)+'\t chi cuadrado ='+str(chicuadrado_G0)+' SI EXCLUIDO\n')
        else:
            archivog.write('g='+str(g)+'\t m='+str(m)+'\t chi cuadrado ='+str(chicuadrado_G0)+' NO EXCLUIDO\n')
            print('NO EXCLUIDO')

archivog.close()
# for masa in masas:
#     m=masa
#     for cte in gtildes:
#         gtilde=cte
#         print()
#         for i in range(nbins):
#             cuentas_conF0[i]=sci.quad(powerlawconalps_F0,bineado[i],bineado[i+1])[0]
#             chicuadrado_F0=np.sum((cuentas_sin-cuentas_conF0)**2/errores**2)
#         print('SIGUIENTE MODELO')
#         print('gtilde=',gtilde)
#         print('masa =', m)
#         print('chi cuadrado F0 = ',chicuadrado_F0)
#         print('treshold = ',treshold)
#         if chicuadrado_F0>=treshold:
#             print('Excluido')
#             archivof.write(str(gtilde)+','+str(m)+','+str(chicuadrado_F0)+',1\n')
#         else:
#             archivof.write(str(gtilde)+','+str(m)+','+str(chicuadrado_F0)+',0\n')
#             print('NO EXCLUIDO')

archivof.close()
# B=5
# errores=cuentas_sin*0.1
# rutaf=r'C:\Users\Usuario\Desktop\pthones\RESULT_F0final_B5.txt'
# archivof=open(rutaf,'a')

# for masa in masas:
#     m=masa
#     for cte in gtildes:
#         gtilde=cte
#         print()
#         for i in range(nbins):
#             cuentas_conF0[i]=sci.quad(powerlawconalps_F0,bineado[i],bineado[i+1])[0]
#             chicuadrado_F0=np.sum((cuentas_sin-cuentas_conF0)**2/errores**2)
#         print('SIGUIENTE MODELO')
#         print('gtilde=',gtilde)
#         print('masa =', m)
#         print('chi cuadrado F0 = ',chicuadrado_F0)
#         print('treshold = ',treshold)
#         if chicuadrado_F0>=treshold:
#             print('Excluido')
#             archivof.write(str(gtilde)+','+str(m)+','+str(chicuadrado_F0)+',1\n')
#         else:
#             archivof.write(str(gtilde)+','+str(m)+','+str(chicuadrado_F0)+',0\n')
#             print('NO EXCLUIDO')

# archivof.close()

# B=25
# # errores=cuentas_sin*0.1
# rutaf=r'C:\Users\Usuario\Desktop\pthones\RESULT_F0final_B25.txt'
# archivof=open(rutaf,'a')

# for masa in masas:
#     m=masa
#     for cte in gtildes:
#         gtilde=cte
#         print()
#         for i in range(nbins):
#             cuentas_conF0[i]=sci.quad(powerlawconalps_F0,bineado[i],bineado[i+1])[0]
#             chicuadrado_F0=np.sum((cuentas_sin-cuentas_conF0)**2/errores**2)
#         print('SIGUIENTE MODELO')
#         print('gtilde=',gtilde)
#         print('masa =', m)
#         print('chi cuadrado F0 = ',chicuadrado_F0)
#         print('treshold = ',treshold)
#         if chicuadrado_F0>=treshold:
#             print('Excluido')
#             archivof.write(str(gtilde)+','+str(m)+','+str(chicuadrado_F0)+',1\n')
#         else:
#             archivof.write(str(gtilde)+','+str(m)+','+str(chicuadrado_F0)+',0\n')
#             print('NO EXCLUIDO')

# archivof.close()

# for masa in masas:
#     m=masa
#     for cte in g1s:
#         g1=cte
#         for cte2 in g2s:
#             g2=cte2
#             for i in range(nbins):
#                 cuentas_conG02[i]=sci.quad(powerlawconalps_G0cuadrado,bineado[i],bineado[i+1])[0]
#                 chicuadrado_G02=np.sum((cuentas_sin-cuentas_conG02)**2/errores**2)
#             print('SIGUIENTE MODELO')
#             print('g1=',g1)
#             print('g2=',g2)
#             print('masa =', m)
#             print('chi cuadrado G02 = ',chicuadrado_G02)
#             print('treshold = ',treshold)
#             if chicuadrado_G02>=treshold:
#                 print('Excluido')
#                 archivog2.write(str(g1)+','+str(g2)+','+str(m)+','+str(chicuadrado_G02)+',1\n')
#             else:
#                 archivog2.write(str(g1)+','+str(g2)+','+str(m)+','+str(chicuadrado_G02)+',0\n')
#                 print('NO EXCLUIDO')



archivog2.close()


print("tiempo de ejecucion")
print("--- %s seconds ---" % (time.time() - start_time))
# print('SIGUIENTE MODELO')
# print('g=',g)
# print('g tilde =',gtilde)
# print('g1=',g1)
# print('g2=',g2)
# print('chi cuadrado G0 = ',chicuadrado_G0)
# print('chi cuadrado F0 = ',chicuadrado_F0)
# print('chi cuadrado G0 cuadrado = ',chicuadrado_G02)
# print('treshold = ',treshold)
# print('sigma = ',stats.chi2.mean(nbins-1))
#####Se podria hacer para que printee automaticamente si se excluye o no, por ahora queda asi y se compara

#####PLOT DEL FLUJO (no necesario)

# plt.figure(figsize=(16,9))
# #no me gusta mucho como quedan las barras de error, si tal cambiar por esto
# plt.plot(bineado[1:],bineado[1:]**2*cuentas_sin,label='No ALPS')
# # plt.errorbar(bineado[1:],bineado[1:]**2*cuentas_sin,yerr=bineado[1:]**2*cuentas_sin*0.05,ecolor='black',label='sin ALPS')
# #ploteo el resto
# plt.plot(bineado[1:],bineado[1:]**2*cuentas_conG0,label=u'w/ALPs, $G_0$ model')
# plt.plot(bineado[1:],bineado[1:]**2*cuentas_conF0,label=u'w/ALPs, $F_0$ model')
# plt.plot(bineado[1:],bineado[1:]**2*cuentas_conG02,label=u'w/ALPs, $G_0^2$ model')
# plt.xscale('log')
# plt.yscale('log')
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.ylabel(u'$E^2\,dN/dE (TeVcm^{-2}s^{-1}$)',fontsize=20)
# plt.xlabel('E(TeV)',fontsize=20)
# plt.legend(fontsize=20)


