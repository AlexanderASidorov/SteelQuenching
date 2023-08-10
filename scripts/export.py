#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 10:50:46 2023

@author: aas
"""

from alloy import Alloy
from factors import Factors

     
def text ():
    Material=Alloy()
    R=Factors.R
    Q=Factors.Q
    fh = open ('data.txt', 'w')
    fh.write('//********************************************************** \n')
    fh.write('//***Coefficients necessary to define model in FlowVision*** \n')
    fh.write('//********************************************************** \n')
    fh.write('\n')
    fh.write('gs=%d; // Grain size coefficient \n' %Material.gs)
    fh.write('R=%.6f;   // Gas constant\n' %R)
    fh.write('Q=%.6f;   // Activation energy \n' %Q)
    fh.write('\n')
    fh.write('//Temperatures of transformation start: \n')
    fh.write('Ae3=%.2f;  // Ferrire transformation start temperature \n' % Material.Ae3)
    fh.write('Ae1=%.2f;  // Pearlite transformation start temperature \n' % Material.Ae1)
    fh.write('Bs=%.2f;   // Bainite transformation start temperature \n' % Material.Bs)
    fh.write('Ms=%.2f;   // Martensite transformation start temperature \n' % Material.Ms)
    fh.write('\n')
    fh.write('//Composition coefficients: \n')
    fh.write('FC=%.7f;  // Ferrire composition coefficient \n' % Material.FC)
    fh.write('PC=%.7f;  // Pearlite composition coefficient \n' % Material.PC)
    fh.write('BC=%.7f;  // Bainite composition coefficient \n' % Material.BC)
    fh.write('MC=%.7f;  // Martensite composition coefficient \n' % Material.alpha_martensite)
    fh.write('\n')
    fh.write('//n coefficients \n')
    fh.write('//Ferrite: \n')
    fh.write('n1_F = %.2f; \n' % Material.n1_F)
    fh.write('n2_F = %.2f; \n' % Material.n2_F)
    fh.write('//Pearlite: \n')
    fh.write('n1_P = %.2f; \n' % Material.n1_P)
    fh.write('n2_P = %.2f; \n' % Material.n2_P)
    fh.write('//Bainite: \n')
    fh.write('n1_B = %.2f; \n' % Material.n1_B)
    fh.write('n2_B = %.2f; \n' % Material.n2_B)
        
    fh.close()
    return fh

def results2text(data):
    fh = open ('results.txt', 'w')
    print(data.iloc[-1, :], file=fh)
    fh.close()
    return fh


  