#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 16:29:15 2023

@author: aas
"""


# Steel 60C2H
Composition={'GrainSize': 7.0, 'C': 0.8, 'Mn': 0.7, 'Si': 1.9, 'Ni': 0.2, 
              'Cr': 0.3, 'Mo': 0.0, 'W': 0.0, 'As': 0.0, 'V': 0, 'Cu': 0.03}
# Test Composition
# Composition={'GrainSize': 5.0, 'C': 0.65, 'Mn': 0.7, 
#         'Si': 1.8, "Ni": 0.2, 'Cr': 0.3, 'Mo': 0.0, 'W': 0.0, 'As':0.0, 'V': 0,
#         'Cu': 0.03}


t, T = [0, 60], [800, 50]

# fraction step (0.01 recomended)
delta_x=0.025
# time step (0.1 - 1)
delta_t=0.1
# Temperature step (1 Celcius recomended)
delta_T=1
# F
# Function to be used (S(x) or I(x))
FunctionType='S(x)'



