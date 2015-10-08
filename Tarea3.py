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