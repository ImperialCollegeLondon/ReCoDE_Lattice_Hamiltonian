# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 16:18:24 2024

@author: ljh3218
"""

import numpy as np
from itertools import product
import scipy.constants as const

def Laguerre(alpha, n, x):
    L_0 = 1
    L_1 = 1 + alpha - x
    if n == 0:
        return L_0
    elif n == 1:
        return L_1
    else:
        L = (1/(n))*((2*n - 1 + alpha - x)*Laguerre(alpha, n-1, x) - (n+alpha-1)*Laguerre(alpha, n-2, x))
        return L

def FCWD_single_nm_v2(n, m, Ereor_inner, E_peak, Ereor_outer, DeltaE):
    #Calculates decay from excited state in vibrational level m to the ground state in vibrational level n
    #Using expression in https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.031055
    kT = 0.0257
    #Huang-Rhys factor
    S = Ereor_inner/E_peak
    #Vibronic integral for any n,m
    prefactor = np.exp(-S)*S**(n-m)*np.math.factorial(m)/np.math.factorial(n)
    lag = Laguerre(n-m, m, S)
    activation = np.exp(-(DeltaE - (n-m)*E_peak - Ereor_outer)**2/(4*Ereor_outer*kT))
    #prefactor is for normalisation - see notebook 3, 14/11/23
    thermal_pop = (1-np.exp(-E_peak/kT))*np.exp(-(m*E_peak)/kT)
    return prefactor*lag**2*activation*thermal_pop

def calc_FCWD_total(Ereor_inner, E_peak, Ereor_outer, DeltaE, N = 20, M = 6):
    kT = 0.0257 
    FCWD_total = 0
    for n, m in product(range(N), range(M), repeat = 1):
        FCWD_total += FCWD_single_nm_v2(n, m, Ereor_inner, E_peak, Ereor_outer, DeltaE)
    #Factor of 1/e to convert eV to joules 
    return (1/const.e)*(1/np.sqrt(4*np.pi*Ereor_outer*kT))*FCWD_total

def decay_rate(Ereor_inner, E_peak, Ereor_outer, DeltaE, V):
    #Convert V into joules
    #Should be a factor of DeltaE as couplings proportional to this
    #But made things tricky for compatiiblity with MLJ only model so have neglected this here
    V *= const.e
    FCWD_0 = calc_FCWD_total(Ereor_inner, E_peak, Ereor_outer, DeltaE)
    return (2*np.pi*V**2*FCWD_0)/const.hbar
    
