#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 17:46:51 2023

@author: aas

In this type of calculatin we use the following modules:
    # solve_JMAK.py - for calculation coefficients of Avrami equation
    # read_data.py - for read TTT data from Excel and for calculation coefficients 
    of Avrami equation via solve_JMAK methods
    # equilibria.py - for calculation of equilibria phase fraction in the case 
    Ferrite and Pearlite reaction can happen simultaneously
    # export.py – for writing necessary for FlowVision software data into text 
    files
    # plot.py - for plotting TTT and plotting example of isothermal 
    transformation
    
"""

from export import tau_n_b_2Excel, tau_2text, n_2text, eq_2text, time_s_2text, phasa_points #, b_2text
from plot import Plot
#from solve_JMAK import JMAK
# from read_data import Create_variables_from_Excel as create
from equilibria import Equilibria


#import os.path
#import pandas as pd
#import numpy as np
#from scipy.optimize import curve_fit
#import math


data=Equilibria('TTT.xlsx', '../data')
Ferrite=data.Ferrite
Pearlite=data.Pearlite
Bainite=data.Bainite


tau_n_b_2Excel(Ferrite, 'Ferrite.xlsx', 'Ferrite')
tau_n_b_2Excel(Pearlite, 'Pearlite.xlsx', 'Pearlite')
tau_n_b_2Excel(Bainite, 'Bainite.xlsx', 'Bainite')

# Export tau, n, start_temperature and equilibrion into a text file in a 
# suitable for FlowVision format
tau_2text(Ferrite, 'Ferrite_tau.txt')
n_2text(Ferrite, 'Ferrite_n.txt')
eq_2text(Ferrite, 'Ferrite_eq.txt')
time_s_2text(Ferrite, 'Ferrite_time_s.txt')

tau_2text(Pearlite, 'Pearlite_tau.txt')
n_2text(Pearlite, 'Pearlite_n.txt')
eq_2text(Pearlite, 'Pearlite_eq.txt')
time_s_2text(Pearlite, 'Pearlite_time_s.txt')

tau_2text(Bainite, 'Bainite_tau.txt')
n_2text(Bainite, 'Bainite_n.txt')
eq_2text(Bainite, 'Bainite_eq.txt')
time_s_2text(Bainite, 'Bainite_time_s.txt')


phasa_points(data)

Plot.plot_f(Pearlite, 600)
Plot.plot_TTT_imp(Ferrite, Pearlite, Bainite)

