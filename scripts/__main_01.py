#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 13:47:50 2023

@author: Alexander Sidorov
Calculate TTT, phase change and phase fractions at defined cooling rate 

"""

from initial_data import Composition, t, T, delta_x, delta_T, delta_t
from solve_equations import Solve
from plot import Plot
from export import text, results2text, TTT2PandasDataframe, TTT2Excel



Steel=Solve() # create an object of the Solve type class
data=Steel.data # write main results into pandas dataframe
TTT=Plot() # plot TTT diagramm
text('data_max_plus_grain_size_7.txt') # write data necessary for FlowVision simulation
results = results2text(data, 'results_max_plus_grain_size_7.txt') # write final phase composition
Plot.plotting_phase_change(data) # plot phase change per time and Temperature
TTT2Excel() # create Excel file with TTT data
TTT_dataframe=TTT2PandasDataframe() # write TTT data into Pandas Dataframe

# tau_f=Steel.Factors.tau_f
# tau_p=Steel.Factors.tau_p
# tau_b=Steel.Factors.tau_b