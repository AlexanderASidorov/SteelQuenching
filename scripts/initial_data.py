#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 16:29:15 2023

@author: aas
"""



Composition={'GrainSize': 7.0, 'C': 0.8, 'Mn': 0.8, 
       'Si': 0.15, "Ni": 0.04, 'Cr': 0.98, 'Mo': 0.21, 'W': 0.0, 'As':0.0, 'V': 0,
       'Cu': 0}

t, T = [0, 120], [800, 20]

# fraction step (0.01 recomended)
delta_x=0.05
# time step (1 second reccomended)
delta_t=0.1
# Temperature step (1 Celcius recomended)
delta_T=1



#alloy = Alloy(gs=7, C=0.4, Mn=0.4, Si=0.15, Ni=0.04, Cr=0.98, Mo=0.21)
