# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 19:29:49 2024

@author: Ayman
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from sympy import symbols, diff, exp, sin, pi, Function

from functions import *
#from resolution_schema_implicite import *

# Paramètres physiques et de discrétisation
Deff = 10**(-10)      # Coefficient de diffusion (m^2/s)
S    = 8 * 10**(-9)   # Constante (mol/m^3/s)
Ce   = 12             # mol/m^3
R    = 0.5            # Longueur spatiale : Rayon
t    = 10**11         # Temps total
Ntot = 1000           # Nombre de points de grille spatiale
Nt   = 100            # Nombre de pas de temps
dr   = R / (Ntot - 1) # Pas d'espace
dt   = t / Nt         # Pas de tempss

#Alternative(critère de liaison entre dt et dr)
#dt = 100000*dr

# Résolution du problème avec un schéma centré pour dérivée spatiale ordre 2 et shéma avant pour dérivée spatiale ordre 1
C_avant_ordre1   = resolution_implicite_avant_ordre1(Ntot, Nt, Deff, S, Ce, dr, dt, R)

# Résolution du problème avec un schéma centré pour dérivée spatiale ordre 2 et shéma centré pour dérivée spatiale ordre 1
C_centree_ordre1 = resolution_implicite_centree_ordre2(Ntot, Nt, Deff, S, Ce, dr, dt, R)

# Initialisation de la grille et des vecteurs à utiliser
r = np.linspace(0, R, Ntot)          #Vecteurs spatial

plt.plot(r, C_avant_ordre1)
plt.xlabel('Position radiale (r)')
plt.ylabel('Concentration (C)')
plt.title('Évolution temporelle de la concentration (B)')
plt.show()

# Tracé du résultat numérique
plt.plot(r, C_centree_ordre1)
plt.xlabel('Position radiale (r)')
plt.ylabel('Concentration (C)')
plt.title('Évolution temporelle de la concentration (F)')
plt.show()

# Tracé du résultat analytique
    ###Evalaution de la solution analytique
C_analytique = solution_analytique(Ntot, r, R, S, Deff, Ce)
    ###plot
plt.plot(r, C_analytique)
plt.xlabel('Position radiale (r)')
plt.ylabel('Concentration (C)')
plt.title('Évaluation analytique de la concentration')
plt.show()

# Calcul des erreurs
    # L1
L1_avant_ordre1     = error_L1(C_avant_ordre1, C_analytique, dr, R)
L1_centree_ordre1   = error_L1(C_centree_ordre1, C_analytique, dr, R)
    # L2
L2_avant_ordre1     = error_L2(C_avant_ordre1, C_analytique, dr, R)
L2_centree_ordre1   = error_L2(C_centree_ordre1, C_analytique, dr, R)
    # Linf
Linf_avant_ordre1   = error_Linf(C_avant_ordre1, C_analytique)
Linf_centree_ordre1 = error_Linf(C_centree_ordre1, C_analytique)

############################################################################
# Etude par rafinement (E)
Ntot_raf = [100, 200, 500, 700]
    # Calcul du vecteur ordre de convergence selon le nombre de rafinements effectues
ordre_convergence_avant = refinement_study_order(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[0]
    # Trace des erreurs
L1     = refinement_study_order(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[1]  
L2     = refinement_study_order(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[2] 
L_inf  = refinement_study_order(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[3] 
pas_dr = refinement_study_order(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[4]


plt.loglog(pas_dr, L1, label='L1')
plt.loglog(pas_dr, L2, label='L2')
plt.loglog(pas_dr, L_inf, label='L_inf')

plt.xlabel("Pas d'espace (dr)")
plt.ylabel('Erreur (C)')
plt.title("Évaluation des erreurs selon le pas d'espace (Étude de raffinement question E)")
plt.legend(["L1", "L2", "Linf"])  # Ajoute une légende avec les labels spécifiés ci-dessus
plt.show()


############################################################################
pause()
############################################################################

# Etude par rafinement (F)
Ntot_raf = [100, 200, 500, 700]
    # Calcul du vecteur ordre de convergence selon le nombre de rafinements effectues
ordre_convergence_center = refinement_study_order_center(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[0]
    # Trace des erreurs
L1     = refinement_study_order_center(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[1]  
L2     = refinement_study_order_center(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[2] 
L_inf  = refinement_study_order_center(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[3] 
pas_dr = refinement_study_order_center(Ntot_raf, Nt, r, R, dr, dt, S, Deff, Ce)[4]

plt.loglog(pas_dr, L1, label='L1')
plt.loglog(pas_dr, L2, label='L2')
plt.loglog(pas_dr, L_inf, label='L_inf')

plt.xlabel("Pas d'espace (dr)")
plt.ylabel('Erreur (C)')
plt.title("Évaluation des erreurs selon le pas d'espace (Étude de raffinement question F)")
plt.legend(["L1", "L2", "Linf"])  # Ajoute une légende avec les labels spécifiés ci-dessus
plt.show()
    