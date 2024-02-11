"""
MEC8211 - Devoir 1 : Verification de code
Fichier : devoir1_post_results.py
Description : Fichier secondaire 
              (a utiliser conjointement avec devoir1_main.py)
Auteur.e.s :
Date de creation du fichier : 10 février 2024
"""

#%% Importation des modules
import numpy as np
import matplotlib.pyplot as plt
import os


""" Post Results """

# Solution

def plot_stationnary_compar(r_l, st_sol, sim_sol,
                            plotting = 'False',
                            path_save = '',
                            title = ''):
    """ Plot the spatial distribution of salt concentration of 
    the stationnary solution and the finite differences solution
    Entrees:
        - r_l:ARRAY of Spatial Nodes
        - st_sol: ARRAY of Stationnary Solution Values at Nodes
        - sim_sol: ARRAY of Simulation Solution Values at Nodes
        - plotting: BOOL to Determine if We Want to Plot Here the Graph
        - path_save: STR Desired Directory When Saving The Graph
        - title: STR Desired Title of the Graph When Saving The File    
    
    Sortie:
        - FIGURE Graphique de la concentration en sel dans le cylindre
        """
        
    plt.figure()
    plt.plot(r_l, st_sol, label = 'Solution Analytique')
    plt.plot(r_l, sim_sol[-1], label = 'Solution par Différences Finies')
    plt.xlabel("Rayon du Cylindre (m)")
    plt.ylabel("Concentration en sel (mol/m^3)")
    plt.title("Comparaison des solutions analytique et numérique")
    plt.legend()
    plt.grid(linestyle = '--')
    if path_save != '':
        os.chdir(path_save)
        if title != '':
            plt.savefig(title+".png", dpi=600)
        else:
            plt.savefig("ComparaisonSolAnalytique.png", dpi=600)
            
    if plotting == True:
        plt.show()


# Erreurs

def convergence_compar(norm_l, n_l, 
                       typAnalyse = "Spatial", 
                       path_save = '',
                       title = ''):

    """ Construit et affiche un graphe de convergences des erreurs selon
    les différentes normes utilisées et enregistre le graphique dans un
    dossier spécifié
    Entrees:
        - norm_l:ARRAY of Tuples with the norm's name'
        - n_l: ARRAY of Abscisses for the Convergence Graph
        - typAnalyse: STR Type of Convergence Analysis
        - path_save: STR Desired Directory When Saving The Graph
        - title: STR Desired Title of the Graph When Saving The File    
    
    Sortie:
        - FIGURE Graphique de la concentration en sel dans le cylindre
        """
        
    plt.figure()
    for name_norm, norm_l in norm_l:
        
        plt.loglog(n_l, norm_l, "s-", label = name_norm)
        
    if typAnalyse == "Spatial":
        plt.title("Convergence Spatiale: Évolution des Erreurs dans le cylindre")
        plt.xlabel(r"$\Delta r$ [m]")
        
    if typAnalyse == "Temporal":    
        plt.title("Convergence Temporelle: Évolution des Erreurs dans le cylindre")
        plt.xlabel(r"$\Delta t$ [s]")
        
    plt.ylabel(r"Erreur [mol/m$^3$]")
    plt.grid(linestyle = '-')
    plt.legend()
    if path_save != '':
        os.chdir(path_save)
        if title != '':
            plt.savefig(title+".png", dpi=600)
        else:
            plt.savefig(typAnalyse + "_Convergence.png", dpi=600)
    plt.show()