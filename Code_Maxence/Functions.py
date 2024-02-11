# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 16:49:57 2024

@author: Maxence
"""
import numpy as np


""" Auxiliary Functions """

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

""" Errors Definition """

def L1(vect):
    
    return np.linalg.norm(vect, ord = 1)

def L2(vect):
    
    return np.linalg.norm(vect, ord = 2)

def L_inf(vect):
    
    return np.linalg.norm(vect, np.inf)