# MEC8211_Hiver2024_Devoir1

This is the repository for the 1st MEC8211 homework of the Winter 2024 term.

## Authors 

Amishga Alphonius, Ayman Benkiran, and Maxence Farin.

# Project Description

This codes solves the stationary state salt concentration (C) in an under-water cylindrical concrete pillar. The problem is modeled as follows:

dC/dT = D_{eff} \nabla^2 C + S

where D_{eff} is the effective diffusion coefficient of salt and S is the source term from a reaction between the salt and a concrete component.

Assuming that the pillar is tall and that the environment is homogenous, through polar coordinates, the problem can be written as a 1D problem along the r-axis.

# Code requirement
This project uses a Python code to solve the problem described above using the finite difference method.

A simple call of the main script (`python3 devoir1_main.py`) produces a convergence analysis of both a first-order spatial finite difference method and a second-order one.

To run the code, the following modules are required:

- matplotlib
- numpy
- scipy 
- os

# Code Architecture

- All source files are located in the `src` folder.
  
  - `devoir1_main.py`: Main code. Can be called using `python3 devoir1_main.py`.
  - `devoir1_functions.py`: Function library for the code.
  - `devoir1_postresults.py`: Postprocessing library for the code.
  - `devoir1_tests_unitaires.py`: Unit test script. Can be called using `python3 devoir1_tests_unitaires.py`
  
- Produced raw solutions are stored in `CSV` format in the `data` folder.
- Postprocessed results (figures and data) are stored in the `results` folder.
- A change log is available in the `doc` folder.

# Code limitations
- At the moment, the code is only written to model a reaction with a constant source term (0th-order reaction). Eventually the code will be able to model a 1st-order reaction.
- At the moment, only the stationary state problem has been verified.





