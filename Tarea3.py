# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:29:14 2015

Resolución de EDOs utilizando Runge Kutta
"""

import math
import numpy as np
import matplotlib.pyplot as plt

u=1.534

def funcion_osc_Van_Der_Pol(y,z): #y'=z
    '''
    Oscilador de Van der Pol reducido a dos EDOs
    '''
    return z, (-y) - (u* ((y**2)-1)*z)
    
def k1(y_n,z_n,h,f):
    f_evaluado = f(y_n,z_n)
    return h * f_evaluado[0], h * f_evaluado[1]  
    
def k2(y_n,z_n,h,f):
    k_1=k1(y_n,z_n,h,f)
    f_evaluado = f(y_n+(h/2.),z_n+(k_1[1]/2.))
    return h * f_evaluado[0], h * f_evaluado[1]
    
def k3(y_n,z_n,h,f):
    k_1=k1(y_n,z_n,h,f)
    k_2=k2(y_n,z_n,h,f)
    f_evaluado = f(y_n+(h),z_n-k_1[1]+2*(k_2[1]))
    return h * f_evaluado[0], h * f_evaluado[1]   

def Runge_Kutta_3(y_n,z_n,h,f):
    '''
    Runge Kutta de 3er orden
    Utiliza los parámetros k1,k2,k3
    '''
    k_1=k1(y_n,z_n,h,f)
    k_2=k2(y_n,z_n,h,f)
    k_3=k3(y_n,z_n,h,f)
    y_n1=y_n+(1./6.)*(k_1[0]+(4.*k_2[0])+k_3[0])
    z_n1=z_n+(1./6.)*(k_1[1]+(4.*k_2[1])+k_3[1])
    return y_n1,z_n1
    
N_steps = 1000
max=20.*math.pi
h = max / N_steps
y = np.zeros(N_steps)
z = np.zeros(N_steps)

y[0] = 0.1
z[0] = 0

i=1
while i<N_steps:
    y[i], z[i] = Runge_Kutta_3(y[i-1], z[i-1], h, funcion_osc_Van_Der_Pol)    
    i+=1
t_rk = [h * j for j in range(N_steps)]

plt.figure(1)
plt.plot(t_rk, z)
plt.xlabel('s',fontsize=18)
plt.ylabel('y', fontsize=18)
plt.figure(2)
plt.plot(y, z)
plt.xlabel('y', fontsize=18)
plt.ylabel('dy/ds', fontsize=18)

y = np.zeros(N_steps)
z = np.zeros(N_steps)

y[0] = 4.0 #nuevamente, esta vez para y0=4
z[0] = 0

i=1
while i<N_steps:
    y[i], z[i] = Runge_Kutta_3(y[i-1], z[i-1], h, funcion_osc_Van_Der_Pol)    
    i+=1

plt.figure(3)
plt.plot(t_rk, z)
plt.xlabel('s',fontsize=18)
plt.ylabel('y', fontsize=18)
plt.figure(4)
plt.plot(y, z)
plt.xlabel('y', fontsize=18)
plt.ylabel('dy/ds', fontsize=18)

sigma=10.
beta=8./3.
ro=28.

def funcion_atractor_lorenz(f,s):
    '''
    Función que recibe un arreglo como parámetro en el cual
    f[0]: x
    f[1]: y
    f[2]: z
    f[3]: dx/ds
    f[4]: dy/ds
    f[5]: dz/ds
    '''
    dfdt=f
    dfdt[0]=f[3]
    dfdt[1]=f[4]
    dfdt[2]=f[5]
    dfdt[3]=sigma*(f[1]-f[0])
    dfdt[4]=f[0]*(ro-f[2])-f[1]
    dfdt[5]=f[0]*f[1]-beta*f[2]
    return dfdt