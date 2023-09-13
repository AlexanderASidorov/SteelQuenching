#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 17:46:51 2023

@author: aas
"""

from plot import Plot
from solve_JMAK import JMAK
from read_data import Create_variables_from_Excel as create
from equilibria import Equilibria
from export import k_n_2Excel, k_2text, n_2text, eq_2text

import os.path
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import math



data=Equilibria('TTT.xlsx', '../data')
Ferrite=data.Ferrite
Pearlite=data.Pearlite
Bainite=data.Bainite


k_n_2Excel(Ferrite, 'Ferrite.xlsx', 'Ferrite')
k_n_2Excel(Pearlite, 'Pearlite.xlsx', 'Pearlite')
k_n_2Excel(Bainite, 'Bainite.xlsx', 'Bainite')

# Export k, n and equilibrion into a text file in a suitable for 
# FlowVision format
k_2text(Ferrite, 'Ferrite_k.txt')
n_2text(Ferrite, 'Ferrite_n.txt')
eq_2text(Ferrite, 'Ferrite_eq.txt')

k_2text(Pearlite, 'Pearlite_k.txt')
n_2text(Pearlite, 'Pearlite_n.txt')
eq_2text(Pearlite, 'Pearlite_eq.txt')

k_2text(Bainite, 'Bainite_k.txt')
n_2text(Bainite, 'Bainite_n.txt')
eq_2text(Bainite, 'Bainite_eq.txt')


#Plot.plot_f(Ferrite, 750)
Plot.plot_TTT_imp(Ferrite, Pearlite, Bainite)

