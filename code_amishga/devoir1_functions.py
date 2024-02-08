"""
MEC8211 - Devoir 1 : Verification de code
Fichier : devoir1_functions.py
Description : Fichier secondaire contenant les fonctions pour le devoir 1
              (a utiliser conjointement avec devoir1_main.py)
Auteur.es :
Date de creation du fichier : 5 février 2024
"""

#%% Importation des modules
import numpy as np

#%% mdf1_rxn_0
"""
Fonction qui resout par le probleme transitoire jusqu'a l'atteinte du regime 
permanent par la methode des differences finies (Schemas d'ordre globaux 1 en 
temps et en espace).
    - En r = 0 : Un schema Gear avant est utilise pour approximer le gradient 
                 de concentration (ordre 2)
    - Pour les points centraux : 
        - Derivee premeire : differentiation avant (ordre 1)
        - Derivee seconde : differentiation centree (ordre 2)
    - En r = R : Une condition de Dirichlet est imposee

Entrees : 
    - prm_prob : Objet qui contient les parametres du probleme
        - C0 : Concentrations initiales [mol/m^3]
        - Ce : Concentration de sel de l'eau salee [mol/m^3]
        - Deff : Coefficient de diffusion effectif de sel dans le beton [m^2/s]
        - R : Rayon du pilier cylindrique [m^3]
        - S : Terme source constant (reaction d'ordre 0) [mol/m^3/s]
    - prm_sim : Objet qui contient les parametres de simulation
        
Sorties : 
    - tf : Temps de fin de simulation [s]
    - c : Vecteur des concentrations au regime permanent [mol/m^3]
"""
def mdf1_rxn_0(prm_prob, prm_sim):
    tf = 0
    diff = 1
    n = prm_sim.n_noeuds
    A = np.zeros((n, n))
    B = np.zeros(n)

    # Condition initiale
    c = np.full(n, prm_prob.C0)
    c[-1] = prm_prob.Ce

    while diff > prm_sim.tol:
        sum_c_prec = sum(c)

        # Gear avant en r = 0
        A[0][0] = 3
        A[0][1] = -4
        A[0][2] = 1
        B[0] = 0

        # Points centraux
        cst1 = prm_sim.dt*prm_prob.Deff
        for i in range(1, n-1):
            cst2 = prm_sim.dr**2 * prm_sim.mesh[i]  # r_i * dr^2
            A[i][i-1] = -cst1*prm_sim.mesh[i]
            A[i][i] = (cst2 + cst1*(prm_sim.dr + 2*prm_sim.mesh[i]))
            A[i][i+1] = -cst1*(prm_sim.dr + prm_sim.mesh[i])
            B[i] = cst2*(c[i] - prm_sim.dt*prm_prob.S)

        # Dirichlet en r = R
        A[-1][-1] = 1
        B[-1] = prm_prob.Ce

        # Resolution du systeme lineraire
        c = np.linalg.solve(A, B)
        tf += prm_sim.dt
        diff = abs(sum(c)-sum_c_prec)/abs(sum_c_prec)
    return tf, c


#%% mdf2_rxn_0
"""
Fonction qui resout par le probleme transitoire jusqu'a l'atteinte du regime 
permanent par la methode des differences finies (Schemas d'ordre globaux 1 en 
temps et 2 en espace).
    - En r = 0 : Un schema Gear avant est utilise pour approximer le gradient 
                 de concentration (ordre 2)
    - Pour les points centraux : 
        - Derivee premiere : differentiation centree (ordre 2)
        - Derivee seconde : differentiation centree (ordre 3)
    - En r = R : Une condition de Dirichlet est imposee

Entrees : 
    - prm_prob : Objet qui contient les parametres du probleme
        - C0 : Concentrations initiales [mol/m^3]
        - Ce : Concentration de sel de l'eau salee [mol/m^3]
        - Deff : Coefficient de diffusion effectif de sel dans le beton [m^2/s]
        - R : Rayon du pilier cylindrique [m^3]
        - S : Terme source constant (reaction d'ordre 0) [mol/m^3/s]
    - prm_sim : Objet qui contient les parametres de simulation
        
Sorties : 
    - tf : Temps de fin de simulation [s]
    - c : Vecteur des concentrations au regime permanent [mol/m^3]
"""
def mdf2_rxn_0(prm_prob, prm_sim):
    tf = 0
    diff = 1
    n = prm_sim.n_noeuds
    A = np.zeros((n, n))
    B = np.zeros(n)

    # Condition initiale
    c = np.full(n, prm_prob.C0)
    c[-1] = prm_prob.Ce

    while diff > prm_sim.tol:
        sum_c_prec = sum(c)

        # Gear avant en r = 0
        A[0][0] = 3
        A[0][1] = -4
        A[0][2] = 1
        B[0] = 0

        # Points centraux
        cst1 = prm_sim.dt*prm_prob.Deff
        for i in range(1, n-1):
            cst2 = 2 * prm_sim.dr**2 * prm_sim.mesh[i]  # 2 * r_i * dr^2
            A[i][i-1] = cst1*(prm_sim.dr - 2*prm_sim.mesh[i])
            A[i][i] = cst2 + 4*cst1*prm_sim.mesh[i]
            A[i][i+1] = -cst1*(prm_sim.dr + 2*prm_sim.mesh[i])
            B[i] = cst2*(c[i] - prm_sim.dt*prm_prob.S)

        # Dirichlet en r = R
        A[-1][-1] = 1
        B[-1] = prm_prob.Ce

        # Resolution du systeme lineraire
        c = np.linalg.solve(A, B)
        tf += prm_sim.dt
        diff = abs(sum(c)-sum_c_prec)/abs(sum_c_prec)
    return tf, c


#%% analytique
"""
Fonction qui calcule la solution analytique aux points du maillage.

Entrees :
    - prm_prob : Objet qui contient les parametres du probleme
        - C0 : Concentrations initiales [mol/m^3]
        - Ce : Concentration de sel de l'eau salee [mol/m^3]
        - Deff : Coefficient de diffusion effectif de sel dans le beton [m^2/s]
        - R : Rayon du pilier cylindrique [m^3]
        - S : Terme source constant (reaction d'ordre 0) [mol/m^3/s]
    - mesh : Vecteur conteant les noeuds (r_i) du probleme 1D [m]
    
Sortie :
    - c : Le profile de concentration radial analytique au régime permanent [mol/m^3]
"""
def analytique(prm_prob, mesh):
    c = [0.25*prm_prob.S/prm_prob.Deff * prm_prob.R**2 * (r**2/prm_prob.R**2 - 1) + prm_prob.Ce for r in mesh]
    return c


#%% erreur_l1
"""
Fonction qui calcule la valeur de l'erreur L1 de la solution numerique obtenue

Entree :
    - c_num = Solution numerique du probleme 1D [mol/m^3]
    - c_analytique = Solution analytique du probleme 1D [mol/m^3]

Sortie :
    - erreur : Norme de l'erreur L1 de la solution numerique [mol/m^3]
"""
def erreur_l1(c_num, c_analytique):
    erreur = sum([abs(ci_num - ci_analytique) for ci_num, ci_analytique in zip(c_num, c_analytique)])
    erreur *= 1/len(c_num)
    return erreur


#%% erreur_l2
"""
Fonction qui calcule la valeur de l'erreur L2 de la solution numerique obtenue

Entree :
    - c_num = Solution numerique du probleme 1D [mol/m^3]
    - c_analytique = Solution analytique du probleme 1D [mol/m^3]

Sortie :
    - erreur : Norme de l'erreur L2 de la solution numerique [mol/m^3]
"""
def erreur_l2(c_num, c_analytique):
    erreur = sum([abs(ci_num - ci_analytique)**2 for ci_num, ci_analytique in zip(c_num, c_analytique)])
    erreur *= 1/len(c_num)
    erreur = np.sqrt(erreur)
    return erreur


#%% erreur_linfty
"""
Fonction qui calcule la valeur de l'erreur L_\infty de la solution numerique obtenue

Entree :
    - c_num = Solution numerique du probleme 1D [mol/m^3]
    - c_analytique = Solution analytique du probleme 1D [mol/m^3]

Sortie :
    - erreur : Norme de l'erreur L_\infty de la solution numerique [mol/m^3]
"""
def erreur_linfty(c_num, c_analytique):
    erreur = max([abs(ci_num - ci_analytique) for ci_num, ci_analytique in zip(c_num, c_analytique)])
    return erreur