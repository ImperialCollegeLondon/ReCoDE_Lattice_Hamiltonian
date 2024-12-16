# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 16:18:24 2024

@author: ljh3218
"""

import numpy as np
from itertools import product
import scipy.constants as const

def laguerre(alpha:int, n:int, x:float) -> float:
    L_0 = 1
    L_1 = 1 + alpha - x
    if n == 0:
        return L_0
    elif n == 1:
        return L_1
    else:
        L = (1/(n))*((2*n - 1 + alpha - x)*laguerre(alpha, n-1, x) - (n+alpha-1)*laguerre(alpha, n-2, x))
        return L

def FCWD_single_nm_v2(n:int, m:int, lambda_inner:float, e_peak:float, lambda_outer:float, w:float) -> float:
    #Calculates decay from excited state in vibrational level m to the ground state in vibrational level n
    #Using expression in https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.031055
    kT = 0.0257
    #Huang-Rhys factor
    S = lambda_inner/e_peak
    #Vibronic integral for any n,m
    prefactor = np.exp(-S)*S**(n-m)*np.math.factorial(m)/np.math.factorial(n)
    lag = laguerre(n-m, m, S)
    activation = np.exp(-(w - (n-m)*e_peak - lambda_outer)**2/(4*lambda_outer*kT))
    #prefactor is for normalisation 
    thermal_pop = (1-np.exp(-e_peak/kT))*np.exp(-(m*e_peak)/kT)
    return prefactor*lag**2*activation*thermal_pop

def calc_FCWD_total(lambda_inner:float, e_peak:float, lambda_outer:float, w:float, N = 20, M = 6) -> float:
    kT = 0.0257 
    FCWD_total = 0
    for n, m in product(range(N), range(M), repeat = 1):
        FCWD_total += FCWD_single_nm_v2(n, m, lambda_inner, e_peak, lambda_outer, w)
    #Factor of 1/e to convert eV to joules 
    return (1/const.e)*(1/np.sqrt(4*np.pi*lambda_outer*kT))*FCWD_total

def decay_rate(lambda_inner:float, e_peak:float, lambda_outer:float, w:float, v:float) -> float:
    #Convert V into joules
    #Should be a factor of DeltaE as couplings proportional to this
    #But only small effect on the couplings for the range of w values considered here
    v *= const.e
    FCWD_0 = calc_FCWD_total(lambda_inner, e_peak, lambda_outer, w)
    return (2*np.pi*v**2*FCWD_0)/const.hbar
    
