# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 10:40:16 2024

@author: Ayman
"""
import numpy as np

from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix

#from resolution_schema_implicite import *

###Function of matrix A(order 1) construction
def coefficients_matrix_construction_order1_approximation(Ntot, Deff, dr, dt, r):
    """
    This function constructs matrix A of coefficients for the implicit scheme.
    """
    # Matrice de coefficients pour le schéma implicite
    A = np.zeros((Ntot,Ntot))
        #Elements internes
    for i in range(1, Ntot - 1):
        A[i, i]     = 1 / dt + 1/dr * (Deff / (r[i]) + 2 * Deff / (dr))
        A[i, i - 1] = -Deff / (dr ** 2) 
        A[i, i + 1] = -Deff / (r[i] * dr) - Deff / (dr ** 2)

    ##Conditions frontieres
        # Première ligne
    # A[0, 0] = -1
    # A[0, 1] = 1
    A[0, 0] = -3
    A[0, 1] = 4
    A[0, 2] = -1
    
        # Dernière ligne
    A[-1, -1] = 1

    return A

###Function of matrix A(order2) construction
def coefficients_matrix_construction_order2_approximation(Ntot, Deff, dr, dt, r):
    """
    This function constructs matrix A of coefficients for the implicit scheme.
    """
    # Matrice de coefficients pour le schéma implicite
    A = np.zeros((Ntot,Ntot))
        #Elements internes
    for i in range(1, Ntot - 1):
        A[i, i]     = 1 / dt + 2 * Deff / (dr ** 2)
        A[i, i - 1] = Deff/(2*r[i]*dr)-Deff / (dr ** 2) 
        A[i, i + 1] = -Deff / (2 * r[i] * dr) - Deff / (dr ** 2)

    ##Conditions frontieres
        # Première ligne
    #A[0, 0] = -1
    #A[0, 1] = 1
    A[0, 0] = -3
    A[0, 1] = 4
    A[0, 2] = -1
    
        # Dernière ligne
    A[-1, -1] = 1
    
    return A

###Function of error L1
def error_L1(numerical_solution, analytical_solution, dx, L):
    """
    This function evaLuates error L1 in 1D
    """
    error = 1/L * np.sum(dx * np.abs(numerical_solution-analytical_solution))
    
    return error

###Function of error L2
def error_L2(numerical_solution, analytical_solution, dx, L):
    """
    This function evaLuates error L2 in 1D
    """
    error = np.sqrt(1/L * np.sum(dx * (numerical_solution-analytical_solution)*(numerical_solution-analytical_solution)))
    
    return error

###Function of error Linf
def error_Linf(numerical_solution, analytical_solution):
    """
    This function evaLuates error L1 in 1D
    """
    error = np.max(np.abs(numerical_solution - analytical_solution))
    
    return error

  
###Function of order of convergence  
def order_convergence(error_coarse, error_refine, h_coarse, h_refine):
    """
    This function evaLuates order of convergence
    """
    order_convergence = np.log(error_coarse/error_refine)/np.log(h_coarse/h_refine) 
    
    return order_convergence

###Function of analytical solution
def solution_analytique(Ntot, r, R, S, Deff, Ce):
    """
    This function evaLuates analytical solution (stationary)
        *R : vector of position
    """
    r = np.linspace(0, R, Ntot)          #Vecteurs spatial
    C_analytique = 0.25 * S/Deff * R**2 * (r**2/R**2 - 1) + Ce
    
    return C_analytique

###Function of refinement study(avant)
def refinement_study_order(vecteur_Ntot, Nt, r, R, dr, dt, S, Deff, Ce):
    """
    This function evaluates the order of convergence using multiple refinement meshes (1D)
    for concentration problem with upwind approximation for order 1
    """
    
    # Initialisation des vecteurs ou scalaires
    Vecteur_order = np.zeros(len(vecteur_Ntot) - 1)
    L1    = np.zeros(len(vecteur_Ntot))
    L2    = np.zeros(len(vecteur_Ntot))
    L_inf = np.zeros(len(vecteur_Ntot))
    vecteur_dr = np.zeros(len(vecteur_Ntot))
    L1_avant_ordre1_raf_temp   = 0
    L2_avant_ordre1_raf_temp   = 0
    Linf_avant_ordre1_raf_temp = 0
    dr_raf_temp                = 0
    
    for i in range(len(vecteur_Ntot)):
        dr_raf   = R / (vecteur_Ntot[i] - 1) # Pas d'espace
            # Resolution
        #dt = 
        C_avant_ordre1_raf   = resolution_implicite_avant_ordre1(vecteur_Ntot[i], Nt, Deff, S, Ce, dr, dt, R) #Solution numerique
        C_analytique_raf     = solution_analytique(vecteur_Ntot[i], r, R, S, Deff, Ce) #Solution analytique
            # Calcul des erreurs
        L1_avant_ordre1_raf     = error_L1(C_avant_ordre1_raf, C_analytique_raf, dr_raf, R)
        L2_avant_ordre1_raf     = error_L2(C_avant_ordre1_raf, C_analytique_raf, dr_raf, R)
        Linf_avant_ordre1_raf   = error_Linf(C_avant_ordre1_raf, C_analytique_raf)
            # Mise a jour des vecteurs des erreurs et pas d'espace
        L1[i]    = L1_avant_ordre1_raf
        L2[i]    = L2_avant_ordre1_raf 
        L_inf[i] = Linf_avant_ordre1_raf
        vecteur_dr[i] = dr_raf
        
        if (i>0):
            # Calcul de l'ordre entre raf_i et raf_(i+1) (en L2)
            Vecteur_order[i-1] = order_convergence(L2_avant_ordre1_raf_temp, L2_avant_ordre1_raf, dr_raf_temp, dr_raf)
    
        # Passation des valeurs pour l'iteration subsequente
        L1_avant_ordre1_raf_temp   = L1_avant_ordre1_raf 
        L2_avant_ordre1_raf_temp   = L2_avant_ordre1_raf
        Linf_avant_ordre1_raf_temp = Linf_avant_ordre1_raf
        dr_raf_temp                = dr_raf
    
    return Vecteur_order, L1, L2, L_inf, vecteur_dr 

###Function of refinement study (centree)
def refinement_study_order_center(vecteur_Ntot, Nt, r, R, dr, dt, S, Deff, Ce):
    """
    This function evaluates the order of convergence using multiple refinement meshes (1D)
    for concentration problem with centered approximation for order 1
    """
    
    # Initialisation des vecteurs ou scalaires
    Vecteur_order = np.zeros(len(vecteur_Ntot) - 1)
    L1    = np.zeros(len(vecteur_Ntot))
    L2    = np.zeros(len(vecteur_Ntot))
    L_inf = np.zeros(len(vecteur_Ntot))
    vecteur_dr = np.zeros(len(vecteur_Ntot))
    L1_centree_ordre1_raf_temp   = 0
    L2_centree_ordre1_raf_temp   = 0
    Linf_centree_ordre1_raf_temp = 0
    dr_raf_temp                  = 0
    
    for i in range(len(vecteur_Ntot)):
        dr_raf   = R / (vecteur_Ntot[i] - 1) # Pas d'espace
            # Resolution    
        C_centree_ordre1_raf   = resolution_implicite_centree_ordre2(vecteur_Ntot[i], Nt, Deff, S, Ce, dr, dt, R) #Solution numerique
        C_analytique_raf     = solution_analytique(vecteur_Ntot[i], r, R, S, Deff, Ce) #Solution analytique

            # Calcul des erreurs
        L1_centree_ordre1_raf     = error_L1(C_centree_ordre1_raf, C_analytique_raf, dr_raf, R)
        L2_centree_ordre1_raf     = error_L2(C_centree_ordre1_raf, C_analytique_raf, dr_raf, R)
        Linf_centree_ordre1_raf   = error_Linf(C_centree_ordre1_raf, C_analytique_raf)
            # Mise a jour des vecteurs des erreurs et pas d'espace
        L1[i]    = L1_centree_ordre1_raf
        L2[i]    = L2_centree_ordre1_raf 
        L_inf[i] = Linf_centree_ordre1_raf
        vecteur_dr[i] = dr_raf
        
        if (i>0):
            # Calcul de l'ordre entre raf_i et raf_(i+1) (en L2)
            Vecteur_order[i-1] = order_convergence(L2_centree_ordre1_raf_temp, L2_centree_ordre1_raf, dr_raf_temp, dr_raf)
    
        # Passation des valeurs pour l'iteration subsequente
        L1_centree_ordre1_raf_temp   = L1_centree_ordre1_raf 
        L2_centree_ordre1_raf_temp   = L2_centree_ordre1_raf
        Linf_centree_ordre1_raf_temp = Linf_centree_ordre1_raf
        dr_raf_temp                  = dr_raf
    
    return Vecteur_order, L1, L2, L_inf, vecteur_dr 

def resolution_implicite_avant_ordre1(Ntot, Nt, Deff, S, Ce, dr, dt, R):
    """
    This function resolves implicite scheme for concentration problem with
    upwind approximation for order 1
    """
    
    # Initialisation de la grille et des vecteurs à utiliser
    r = np.linspace(0, R, Ntot)          #Vecteur spatial
    C = np.zeros(Ntot)                   #Vecteur concentration

    vecteur_S     = np.zeros(Ntot) - S   #Vecteur terme source
    vecteur_S[0]  = 0                
    vecteur_S[-1] = 0

    CI     = np.zeros(Ntot)              #Vecteur conditions initiales
    CI[-1] = Ce

    # Conditions initiales
    C[0]  = 0
    C[-1] = Ce 

    # Construction de la matrice des coefficients du problème
    A = coefficients_matrix_construction_order1_approximation(Ntot, Deff, dr, dt, r)
    A = csc_matrix(A) #Conversion de la matrice A en format csc
    
    # Boucle temporelle pour le schéma implicite
    C_avant = 0
    for it in range(Nt):
        vecteur_B     = (1/dt) * C
        vecteur_B[0]  = 0
        vecteur_B[-1] = 0
            ##Membre de droite
        b = vecteur_B + vecteur_S + CI

            ##Resolution matricielle en temps it
                ###Avec linalg.solve
        #C = np.linalg.solve(A, b)
                ###Avec spsolve
        C = spsolve(A, b)
            ##Critere d'arret
        if np.sum(C_avant) != 0:
            if (abs(np.sum(C_avant - C))/np.sum(C_avant) < 0.0000001):
                break;
        C_avant = C
    
    return C

def resolution_implicite_centree_ordre2(Ntot, Nt, Deff, S, Ce, dr, dt, R):
    """
    This function resolves implicite scheme for concentration problem with
    centered approximation for order 1
    """
    
    # Initialisation de la grille et des vecteurs à utiliser
    r = np.linspace(0, R, Ntot)          #Vecteurs spatial
    C = np.zeros(Ntot)                   #Vecteur concentration

    vecteur_S     = np.zeros(Ntot) - S   #Vecteur terme source
    vecteur_S[0]  = 0                
    vecteur_S[-1] = 0

    CI     = np.zeros(Ntot)              #Vecteur conditions initiales
    CI[-1] = Ce

    # Conditions initiales
    C[0]  = 0
    C[-1] = Ce 

    # Construction de la matrice des coefficients du problème
    A = coefficients_matrix_construction_order2_approximation(Ntot, Deff, dr, dt, r)
    A = csc_matrix(A) #Conversion de la matrice A en format csc
    
    # Boucle temporelle pour le schéma implicite
    for it in range(Nt):
        vecteur_B     = (1/dt) * C
        vecteur_B[0]  = 0
        vecteur_B[-1] = 0
            ##Membre de droite
        b = vecteur_B + vecteur_S + CI
            ##Resolution matricielle en temps it
                ###Avec linalg.solve
        #C = np.linalg.solve(A, b)
                ###Avec spsolve
        C = spsolve(A, b)
        
    return C

def pause():
    """
    This function allows to apply pause of the execution
    """
    
    input('Veuillez cliquer sur entrer pour coninuer')
    