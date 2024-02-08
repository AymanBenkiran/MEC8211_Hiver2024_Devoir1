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
Class qui contient les differents parametres du probleme selon le type de terme source

Entree:
    - ordre_de_rxn : Ordre de la reaction (0 ou 1)
"""
class ParametresProb:
    def __init__(self, ordre_de_rxn=0):
        self.C0 = 0         # Concentrations initiales [mol/m^3]
        self.Ce = 12        # Concentration de sel de l'eau salee [mol/m^3]
        self.R = 1*0.5      # Rayon du pilier cylindrique [m^3]
        self.Deff = 1e-10   # Coefficient de diffusion effectif de sel dans le beton [m^2/s]
        if ordre_de_rxn == 0:
            self.S = 8e-9   # Terme source constant (reaction d'ordre 0) [mol/m^3/s]
        elif ordre_de_rxn == 1:
            self.k = 4e-9   # Constante de réaction pour la reaction d'ordre 1 [s^{-1}]
        else:
            raise Exception("L'ordre de reaction doit etre de 0 ou 1.")

"""
Classe qui contient les parametres de simulation pour une simulation donnee

Entrees:
    - prm_rxn : Object contenant les donnees du probleme
        - R : Rayon du pilier cylindrique [m^3]
    - p_n_noeuds : Nombre de noeuds dans le maillage [noeud]
   
Attributs :
    - n_noeuds : Nombre de noeuds dans le maillage [noeud]
    - dr : Pas en espace des differents maillages [m]
    - dt : Pas de temps des differents maillages [s]
    - mesh : Vecteur conteant les noeuds (r_i) du probleme 1D [m]
    - tol : Tolerance relative pour l'atteinte du regime permanent []
    - c : 
"""
class ParametresSim:
    def __init__(self, prm_rxn, p_n_noeuds):
        self.n_noeuds = p_n_noeuds                          # Nombre de noeuds dans le maillage [noeud]
        self.dr = prm_rxn.R/(self.n_noeuds-1)               # Pas en espace des differents maillages [m]
        self.dt = 0.5*self.dr**2/prm_rxn.Deff               # Pas de temps des differents maillages [s]
        self.mesh = np.linspace(0,prm_rxn.R,self.n_noeuds)  # Vecteur conteant les noeuds (r_i) du probleme 1D [m]
        self.tol = 1e-8                                     # Tolerance relative pour l'atteinte du regime permanent []
        self.c = np.zeros(self.n_noeuds)                    # Solution une fois l'atteinte du regime permanent [mol/m^3]
        self.tf = 0                                         # Temps de fin de la simulation [s]


#%% Initialisation des objects contenants les parametres du probleme pour une reaction d'ordre 0 et 1

prm_rxn_0 = ParametresProb(ordre_de_rxn=0)
prm_rxn_1 = ParametresProb(ordre_de_rxn=1)


#%% Discretisation pour la reaction d'ordre 0

n_noeuds = [10, 20, 40, 80, 160] # Liste de nombre de noeuds pour les differents maillages [noeud]

# Initialisation des differents maillages a l'etude
prm_simulations_mdf1_rxn0 = []
prm_simulations_mdf2_rxn0 = []
for i in range(len(n_noeuds)):
    prm_simulations_mdf1_rxn0.append(ParametresSim(prm_rxn_0, n_noeuds[i]))
    prm_simulations_mdf2_rxn0.append(ParametresSim(prm_rxn_0, n_noeuds[i]))


#%% Resolution du probleme

# for prm_sim in prm_simulations_mdf1_rxn0:
#
#     # Calcul de la concentration au regime permanent
#     tf,prm_sim.c = mdf1_rxn_0(prm_rxn_0, prm_sim)
#
#     # Exportation des solutions dans des fichiers csv
#     exported_data = pd.DataFrame({'r': prm_sim.mesh, 'C(r)': prm_sim.c})
#     exported_data.to_csv(f"./solutions/mdf1_rxn0_noeuds_{str(prm_sim.n_noeuds).zfill(3)}.csv", index=False)
#
#     # Affichage au terminal
#     print(f"****************************************************************************")
#     print(f"Nombre de noeuds : {prm_sim.n_noeuds}")
#     print(f"dr = {prm_sim.dr} m")
#     print(f"dt = {prm_sim.dt} s")
#     print(f"tf = {tf} s")
#     print(f"****************************************************************************")

for prm_sim in prm_simulations_mdf2_rxn0:

    # Calcul de la concentration au regime permanent
    tf,prm_sim.c = mdf2_rxn_0(prm_rxn_0, prm_sim)

    # Exportation des solutions dans des fichiers csv
    exported_data = pd.DataFrame({'r': prm_sim.mesh, 'C(r)': prm_sim.c})
    exported_data.to_csv(f"./solutions/mdf2_rxn0_noeuds_{str(prm_sim.n_noeuds).zfill(3)}.csv", index=False)

    # Affichage au terminal
    print(f"****************************************************************************")
    print(f"Nombre de noeuds : {prm_sim.n_noeuds}")
    print(f"dr = {prm_sim.dr} m")
    print(f"dt = {prm_sim.dt} s")
    print(f"tf = {tf} s")
    print(f"****************************************************************************")


#%% Etude de l'erreur

# Initialisation de listes pour entreposer les donnees pour l'affichage graphique
dr = []
liste_erreur_l1 = []
liste_erreur_l2 = []
liste_erreur_linfty = []

# Calcul des erreurs pour les differentes simualtions
for prm_sim in prm_simulations_mdf2_rxn0:
    dr.append(prm_sim.dr)                                               # Pas en espace [m]
    c_analytique = analytique(prm_rxn_0, prm_sim.mesh)                  # Solution analytique [mol/m^3]
    liste_erreur_l1.append(erreur_l1(prm_sim.c, c_analytique))          # Erreur L1
    liste_erreur_l2.append(erreur_l2(prm_sim.c, c_analytique))          # Erreur L2
    liste_erreur_linfty.append(erreur_linfty(prm_sim.c, c_analytique))  # Erreur L_\infty

# Exportation des valeurs d'erreur dans un fichier csv
exported_data = pd.DataFrame({'dr': dr, 'L1_error': liste_erreur_l1, 'L2_error': liste_erreur_l2, 'Linfty_error': liste_erreur_linfty})
exported_data.to_csv("./erreurs/mdf2_rxn0.csv", index=False)

# Affichage graphique
plt.loglog(dr,liste_erreur_l1, "s-", label=r"Erreur L$^1$")
plt.loglog(dr,liste_erreur_l2, "s-", label=r"Erreur L$^2$")
plt.loglog(dr,liste_erreur_linfty, "s-", label=r"Erreur L$^\infty$")
plt.xlabel(r"$\Delta r$ [m]")
plt.ylabel(r"Erreur [mol/m$^3$]")
plt.title("")
plt.legend()
plt.grid()
plt.savefig("erreurs_mdf2.png", dpi=600)
plt.show()