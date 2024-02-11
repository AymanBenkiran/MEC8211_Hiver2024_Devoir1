# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:27:31 2024

@author: Maxence
"""

import logging
import matplotlib.pyplot as plt
import numpy as np


""" Post Results """

# Solution

def plot_stationnary_compar(r_l, st_sol, sim_sol):
    """ Plot the spatial distribution of salt concentration of 
    the stationnary solution and the finite differences solution
    r_l:Array of Spatial Nodes
    st_sol: Array 
    sim_sol: 
        """
    plt.figure()
    plt.plot(r_l, st_sol, label = 'Solution Analytique')
    plt.plot(r_l, sim_sol[-1], label = 'Solution par Différences Finies')
    plt.xlabel("Rayon du Cylindre (m)")
    plt.ylabel("Concentration en sel (mol/m^3)")
    plt.legend()
    plt.grid(linestyle = '--')
    plt.show()


# Erreurs

def convergence_compar(norm_l, n_l, typAnalyse = "Spatial"):
    
    plt.figure()
    for name_norm, norm_l in norm_l:
        
        plt.loglog(n_l, norm_l, label = name_norm)
        
    if typAnalyse == "Spatial":
        plt.title("Évolution des Erreurs dans le cylindre")
        plt.xlabel("Rayon du Cylindre r (m)")
        
    if typAnalyse == "Temporal":    
        plt.title("Évolution des Erreurs dans le cylindre")
        plt.xlabel("Temps écoulé (s)")
        
    plt.ylabel("Erreurs")
    plt.grid(linestyle = '-')
    plt.legend()
    
    plt.show()