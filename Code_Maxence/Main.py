# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 23:05:22 2024

@author: Maxence
"""

# Bibliothèques pertinentes

import numpy as np
import matplotlib.pyplot as plt
from copy import *
import logging
import os

""" Résolution du problème """

# Definition des constantes

S = 8e-9 # mol/m^3/s
D = 1.0 # m
D_eff = 10e-10 #m^/s
C_e = 12.0 # mol/m^3
k = 4e-9 # s^-1
R = 1.0 # Rayon du cylindre

n = 5 # Nombre de noeuds
tol = 1e-5 # Tolérance sur la solution entre deux étapes
t_sim = 0 # Temps simulé pour atteindre la tolérance requise sur l'approximation

dr = R/(n-1) #Pas spatial uniforme
dt = 100.0 #Pas de temps

import pde_def as pf
import results_processing as respr


# Matrices et vecteurs de calcul

r_i = np.linspace(0.0, R, n)

# Conditions de Dirichlet
C_Dirich = np.array([0.0 for i in range(n)])
C_Dirich[-1] = C_e

im1_a = np.array([-D_eff*dt for i in range(n-1)])

r_i[0] = 1.0
i_a = np.array(coeff_iv(dt, dr, r_i, D_eff))
r_i[0] = 0.0

ip1_a = np.array(coeff_ipv(dt, dr, r_i[1:], D_eff))

A = tridiag(im1_a, i_a, ip1_a)

# Schéma de Gear avant

A[0][0], A[0][1], A[0][2] = -1.5*dt*dr, 2*dt*dr, -0.5*dt*dr

A[-1][-1], A[-1][-2] = dt * dr ** 2, 0

# Solution Stationnaire

def stationnary_sol(S, D_eff, C_e, R, r):
    
    return 0.25*S/D_eff*(r**2-R**2) + C_e

st_sol_v = np.vectorize(stationnary_sol)

st_sol_res = st_sol_v(S, D_eff, C_e, R, r_i) 

#Liste d'erreurs

L1_l = []
L2_l = []
Linf_l = []

# Représentation d'état des concentrations

C0 = np.array([0 for i in range(n)])
C0[-1] = C_e
C = [C0]

C_tm = C[-1]
C_tm[-1] = 0
term_droite = (dr**2) * (C_tm + dt * (S+C_Dirich))
C_t = np.linalg.solve(A, term_droite)
C.append(C_t)
t_sim += dt

L1_l.append(L1(C[-1]-st_sol_res))
L2_l.append(L2(C[-1]-st_sol_res))
Linf_l.append(L_inf(C[-1]-st_sol_res))

while 1e2*L_inf(C[-1]-C[-2]) > tol:
    print(L_inf(C[-1]-C[-2]))
    C_tm = deepcopy(C[-1])
    C_tm[-1] = 0.0
    term_droite = (dr**2) * (C_tm + dt * (S+C_Dirich))
    C_t = np.linalg.solve(A, term_droite)
    C.append(C_t)
    
    t_sim += dt
    L1_l.append(L1(C[-1]-st_sol_res))
    L2_l.append(L2(C[-1]-st_sol_res))
    Linf_l.append(L_inf(C[-1]-st_sol_res))
    
    
""" Analyse des Résultats """
    
# Solution

plt.figure()
plt.plot(r_i, st_sol_res, label = 'Solution Analytique')
plt.plot(r_i, C[-1], label = 'Solution par Différences Finies')
plt.xlabel("Rayon du Cylindre (m)")
plt.ylabel("Concentration en sel (mol/m^3)")
plt.legend()
plt.grid(linestyle = '--')
plt.show()

# Erreurs

n_t = len(L1_l)
t = np.linspace(0.0, t_sim, n_t)

plt.figure()

plt.loglog(t, L1_l, label = 'Norme 1')
plt.loglog(t, L2_l, label = 'Norme 2')
plt.loglog(t, Linf_l, label = 'Norme Infinie')

plt.title("Évolution des Erreurs dans le cylindre")
plt.xlabel("Temps écoulé (s)")
plt.ylabel("Erreurs")
plt.grid(linestyle = '-')
plt.legend()

plt.show()
