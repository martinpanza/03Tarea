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