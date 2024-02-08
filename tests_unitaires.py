# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 14:59:02 2024

@author: Ayman
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from sympy import symbols, diff, exp, sin, pi, Function

from functions import *

def test_errors():
    """
    This function verifies the L1, L2 and L_inf errors between two vectors
    """
    
    #Vectors and scalars of inputs
    x = [1, 2, 3, 4, 5]
    vecteur_x = np.array(x)
    
    y = [1.01, 2.1, 3.2, 4.1, 5.11]
    vecteur_y = np.array(y)
    
    dx = 0.1
    L  = 1.0
    
    #Outputs (execution)
    L1 = error_L1(vecteur_x, vecteur_y, dx, L)
    L2 = error_L2(vecteur_x, vecteur_y, dx, L)
    L_inf = error_Linf(vecteur_x, vecteur_y)
    
    #Verification with expected results
    test1 = 0
    test2 = 0
    test3 = 0
    
        #Conditions
    if (abs(L1 - 0.052) < 0.001) :
        test1 = 1
    if (abs(L2 - 0.08497) < 0.0001):
        test2 = 1
    if (abs(L_inf - 0.2) < 0.0001):
        test3 = 1
    
        #Print test results
    if (test1 == 1):
        print("L'erreur L1 est verifiee")
    elif (test1 == 0):
        print("Il y a une erreur dans l'erreur L1")
        print("La solution trouvee est :")
        print(L1)
        
    if (test2 == 1):
        print("L'erreur L2 est verifiee")

    elif (test2 == 0):
        print("Il y a une erreur dans l'erreur L2")
        print("La solution trouvee est :")
        print(L2)
        
    if (test3 == 1):
        print("L'erreur Linf est verifiee")
    elif (test3 == 0):
        print("Il y a une erreur dans l'erreur Linf")
        print("La solution trouvee est :")
        print(L_inf)
  
def test_order_convergence():
    """
    This function verifies the evaluation of order of convergence
    """
    
    #Inputs
    error_coarse = 0.01
    error_refine = 0.001
    h_coarse = 100
    h_refine = 0.1
    
    #Output (execution)
    ordre = order_convergence(error_coarse, error_refine, h_coarse, h_refine)

    #Verification with expected results
    test = 0
    
        #Conditions
    if (abs(ordre - 0.3333) < 0.0001) :
        test = 1
    
        #Print test results
    if (test == 1):
        print("L'ordre de convergence est verifie")

    elif (test == 0):
        print("Il y a une erreur dans l'ordre de convergence")
        print("La solution trouvee est :")
        print(ordre)
        
def test_solution_analytique():
    """
    This function verifies the evaluation of analytical solution
    """
    
    #Inputs
    Ntot = 100
    R    = 0.5
    r    = np.linspace(0, R, Ntot)
    S    = 8*10**(-9)
    Deff = 10**(-10)
    Ce   = 12
    
    #Output (execution)
    solution = solution_analytique(Ntot, r, R, S, Deff, Ce)

    #Verification with expected results
    test = 0
    
        #Conditions
    if (abs(solution[0] - 7) < 0.0001) and  (abs(solution[-1] - 12) < 0.0001):
        test = 1
    
        #Print test results
    if (test == 1):
        print("La solution analytique en r = 0 et r = R est verifiee")

    elif (test == 0):
        print("Il y a une erreur dans la solution analytique")
        print("La solution trouvee est :")
        print("*en r = 0")
        print(solution[0])
        print("*en r = R")
        print(solution[-1])
        
#Execution de tests unitaires
print("########################")
print("###Verification erreurs")
test_errors()

print("####################################")
print("###Verification ordre de convergence")
test_order_convergence()

print("####################################")
print("###Verification solution analytique")
test_solution_analytique()