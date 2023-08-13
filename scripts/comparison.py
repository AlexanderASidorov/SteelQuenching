#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 13:51:38 2023

@author: aas
"""

from solve_equations import Solve
from plot import Plot
import os.path
import pandas as pd

Steel=Solve() # create an object of the Solve type class
data=Steel.data # write main results into pandas dataframe

os.chdir('../data') # moving to the directory data
FVdata = pd.read_csv('fvresults.txt', sep='\s+')
# os.chdir('../scripts') # moving to the directory scripts

Plot.plotting_comparison(data, FVdata, 'Austenite')
Plot.plotting_comparison(data, FVdata, 'Ferrite')
Plot.plotting_comparison(data, FVdata, 'Pearlite')
Plot.plotting_comparison(data, FVdata, 'Bainite')
Plot.plotting_comparison(data, FVdata, 'Martensite')
