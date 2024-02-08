'''
MEC8211 - Devoir 1 : Verification de code

Fichier : devoir1_main.py
Description : Fichier principal pour le devoir 1
Lancer avec :  python3 devoir1_main.py
Auteur.es :
Date de creation du fichier : 5 fevrier 2024
'''

#%% Importation des modules
import matplotlib.pyplot as plt
import pandas as pd
# import sys args

# Importation des fonctions
try:
    from devoir1_functions import *
except:
    print("ERREUR ! Il y a une erreur fatale dans le fichier devoir1_functions.py")

#%% Donnees du probleme et definition des classes

"""
Classe qui contient les differents parametres du probleme selon le type de terme source

Entree:
    - ordre_de_rxn : Ordre de la reaction (0 ou 1) []
    
Attributs:
    - c0 : float - Concentrations initiales [mol/m^3]
    - ce : float - Concentration de sel de l'eau salee [mol/m^3]
    - r : float - Rayon du pilier cylindrique [m]
    - d_eff : float - Coefficient de diffusion effectif de sel dans le beton [m^2/s]
    - ordre_de_rxn : int - Ordre de la reaction (0 ou 1) []
    - s : float - Terme source constant (reaction d'ordre 0) [mol/m^3/s]
    - k : float - Constante de réaction pour la reaction d'ordre 1 [s^{-1}]
"""
class ParametresProb:
    def __init__(self, ordre_de_rxn=0):
        self.c0 = 0.0
        self.ce = 12.0
        self.r = 1*0.5
        self.d_eff = 1e-10
        self.ordre_de_rxn = ordre_de_rxn
        if self.ordre_de_rxn == 0:
            self.s = 8e-9
        elif self.ordre_de_rxn == 1:
            self.k = 4e-9
            raise Exception("Cet ordre n'est pas encore implememte.")
        else:
            raise Exception("L'ordre de reaction doit etre de 0 ou 1.")

"""
Classe qui contient les parametres de simulation pour une simulation donnee

Entrees:
    - prm_rxn : Object contenant les donnees du probleme
        - r : Rayon du pilier cylindrique [m^3]
    - p_n_noeuds : Nombre de noeuds dans le maillage [noeud]
    - p_mdf : 
   
Attributs :
    - n_noeuds : int - Nombre de noeuds dans le maillage [noeud]
    - dr : float - Pas en espace des differents maillages [m]
    - dt : float - Pas de temps des differents maillages [s]
    - mesh : array of floats - Vecteur conteant les noeuds (r_i) du probleme 1D [m]
    - tol : float - Tolerance relative pour l'atteinte du regime permanent []
    - c : array of floats - Solution une fois l'atteinte du regime permanent [mol/m^3]
    - tf : float - Temps de fin de la simulation [s]
    - mdf : 
    - ordre_de_rxn : 
"""
class ParametresSim:
    def __init__(self, prm_rxn, p_n_noeuds, p_mdf):
        self.n_noeuds = p_n_noeuds
        self.dr = prm_rxn.r/(self.n_noeuds-1)
        self.dt = 0.5*self.dr**2/prm_rxn.d_eff
        self.mesh = np.linspace(0,prm_rxn.r,self.n_noeuds)
        self.tol = 1e-8
        self.c = np.zeros(self.n_noeuds)
        self.tf = 0
        self.mdf = p_mdf
        self.ordre_de_rxn = 0


#%% Initialisation des objects contenants les parametres du probleme pour une reaction d'ordre 0
#   et 1

prm_rxn_0 = ParametresProb(ordre_de_rxn=0)
# prm_rxn_1 = ParametresProb(ordre_de_rxn=1) # pas encore implemente


#%% Discretisation pour la reaction d'ordre 0

n_noeuds = [10, 20, 40, 80, 160] # Liste de nombre de noeuds pour les differents maillages [noeud]

# Initialisation des differents maillages a l'etude
prm_simulations_mdf1_rxn0 = []
prm_simulations_mdf2_rxn0 = []
for i in range(len(n_noeuds)):
    prm_simulations_mdf1_rxn0.append(ParametresSim(prm_rxn_0, n_noeuds[i], 1))
    prm_simulations_mdf2_rxn0.append(ParametresSim(prm_rxn_0, n_noeuds[i], 2))


#%% Resolution du probleme

for prm_simulation in [prm_simulations_mdf1_rxn0, prm_simulations_mdf2_rxn0]:
    for prm_sim in prm_simulation:
        mdf_i = prm_sim.mdf
        # Calcul de la concentration au regime permanent
        exec(f"prm_sim.tf,prm_sim.c = mdf{mdf_i}_rxn_0(prm_rxn_0, prm_sim)")

        # Exportation des solutions dans des fichiers csv
        exported_data = pd.DataFrame({'r': prm_sim.mesh, 'C(r)': prm_sim.c})
        exported_data.to_csv(f"./solutions/mdf{mdf_i}_rxn0_noeuds_"
                             f"{str(prm_sim.n_noeuds).zfill(3)}.csv", index=False)

        # Affichage au terminal
        print("****************************************************************************")
        print(f"Nombre de noeuds : {prm_sim.n_noeuds} noeuds")
        print(f"dr = {prm_sim.dr} m")
        print(f"dt = {prm_sim.dt} s")
        print(f"tf = {prm_sim.tf} s")
        print("****************************************************************************")

#%% Etude de l'erreur

# Initialisation de listes pour entreposer les donnees pour l'affichage graphique
dr = []                     # Pas en espace [m]
liste_erreur_l1 = []        # Erreur L1 [mol/m^3]
liste_erreur_l2 = []        # Erreur L2
liste_erreur_linfty = []    # Erreur Linfty

# Calcul des erreurs pour les differentes simualtions
for prm_simulation in [prm_simulations_mdf1_rxn0, prm_simulations_mdf2_rxn0]:
    dr = []
    liste_erreur_l1 = []
    liste_erreur_l2 = []
    liste_erreur_linfty = []
    for prm_sim in prm_simulation:
        mdf_i = prm_sim.mdf
        dr.append(prm_sim.dr)
        c_analytique = analytique(prm_rxn_0, prm_sim.mesh)  # Solution analytique [mol/m^3]
        liste_erreur_l1.append(erreur_l1(prm_sim.c, c_analytique))
        liste_erreur_l2.append(erreur_l2(prm_sim.c, c_analytique))
        liste_erreur_linfty.append(erreur_linfty(prm_sim.c, c_analytique))

    # Exportation des valeurs d'erreur dans un fichier csv
    exported_data = pd.DataFrame({'dr': dr, 'L1_error': liste_erreur_l1,
                                  'L2_error': liste_erreur_l2,
                                  'Linfty_error': liste_erreur_linfty})
    exported_data.to_csv(f"./erreurs/mdf{mdf_i}_rxn0.csv", index=False)

    # Affichage graphique
    plt.figure()
    plt.loglog(dr,liste_erreur_l1, "s-", label=r"Erreur L$_1$")
    plt.loglog(dr,liste_erreur_l2, "s-", label=r"Erreur L$_2$")
    plt.loglog(dr,liste_erreur_linfty, "s-", label=r"Erreur L$_\infty$")
    plt.xlabel(r"$\Delta r$ [m]")
    plt.ylabel(r"Erreur [mol/m$^3$]")
    plt.title("")
    plt.legend()
    plt.grid()
    plt.savefig(f"erreurs_mdf{mdf_i}.png", dpi=600)
    plt.show()
