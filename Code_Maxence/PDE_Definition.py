# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:20:07 2024

@author: Maxence
"""

import numpy as np
import logging

def coeff_i(dt, dr, r_i, D_eff, order = "1"):
    
    if order == "1":
        
        return dr ** 2 + dt * D_eff * (dr/(r_i) + 2*dt)
    
    elif order == "2":
        
        return dr ** 2 + dt * D_eff * 2 * dt
    else:
        raise Exception("Method Non Implemented for Order Specified")
        
coeff_iv = np.vectorize(coeff_i)

def coeff_ip(dt, dr, r_i, D_eff, order = "1"):
    
    if order == "1":
        return -D_eff * dt * (dr/r_i + 1)
    
    elif order == "2":
        return -D_eff * dt * (0.5 * dr/r_i + 1)
    
    else:
        raise Exception("Method Non Implemented for Order Specified")

coeff_ipv = np.vectorize(coeff_ip)

def coeff_im(dt, dr, r_i, D_eff, order = "1"):
    if order == "1":
        return -D_eff*dt
    
    elif order == "2":
        return -D_eff * dt + 0.5 * D_eff * dt * dr / r_i
    
    else:
        raise Exception("Method Non Implemented for Order Specified")

coeff_imv = np.vectorize(coeff_im)


def terme_source(S, r_i, typ_source = "const"):
    
    n = len(r_i)
    
    if typ_source == "const":
        
        S_vect = np.array([0.0] + [-S for i in range(n-2)] + [0.0])
        
    elif typ_source == "function":
        
        S_l = []
        for r in r_i:
            S_l.append(S(r))
        S_vect = np.array(S_l)
    
    return S_vect