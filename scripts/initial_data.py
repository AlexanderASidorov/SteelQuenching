#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 16:29:15 2023

@author: aas
"""


# Steel 40H (see JMatPro demo)
Composition={'GrainSize': 5.0, 'C': 0.65, 'Mn': 0.5, 
       'Si': 1.45, "Ni": 0.2, 'Cr': 0.25, 'Mo': 0.0, 'W': 0.0, 'As':0.0, 'V': 0,
       'Cu': 0.03}

t, T = [0, 1000], [850, 60]

# fraction step (0.01 recomended)
delta_x=0.025
# time step (0.1 - 1)
delta_t=0.1
# Temperature step (1 Celcius recomended)
delta_T=5
# F
# Function to be used (S(x) or I(x))
FunctionType='S(x)'



