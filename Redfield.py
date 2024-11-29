# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 16:19:11 2024

@author: ljh3218
"""

import numpy as np
import scipy.constants as sp
import scipy.integrate as integrate

#factor of h/e here as we input the reorganisation energies in eV and C(w) should have units of s-1
def C_re_2D_array(w, reorE, Epeak, kT):
    C_w = ((1+nw(w,kT))*((J_Renger(w,reorE,Epeak,kT = kT))\
        - (J_Renger(-w,reorE,Epeak,kT = kT))))/(sp.hbar/sp.e)   
    return 2*np.pi*C_w

def nw(hw,kT):                       
    '''Bose-Einstein occupancy function'''
    if hasattr(hw, "__len__") == 0:
        if hw == 0: n=99999
        else: n=1/(np.exp(hw/kT)-1) 
    elif hasattr(hw, "__len__") == 1:
        hw[hw == 0] = 1e-7
        n=1/(np.exp(hw/kT)-1) 
    return n                                                 

def J_Renger(w, reorE, Epeak, kT = 0.0257, n = 15):
    a = Epeak/(3/n)**(1/n)
    J_w = reorE*w**3*np.exp((-w/a)**n)
    J_w[w < 0] = 0
    return J_w
   
def norm_J_Renger(reorE, Epeak, n = 15):
    def func(w, a, n):
        return w**2*np.exp(-(w/a)**n)
    w_end = 0.3
    a = Epeak/(3/n)**(1/n)
    I0 = integrate.quad(lambda x: func(x, a, n), 0, w_end)
    reorE_input = reorE/I0[0]
    return reorE_input
