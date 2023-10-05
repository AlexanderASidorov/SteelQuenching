#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 10:50:46 2023

@author: aas
"""
import numpy as np
import os.path
import pandas as pd
from alloy import Alloy
from factors import Factors
from fraction import Fraction
from openpyxl import Workbook

from initial_data import Composition, t, T, delta_x, delta_T, delta_t, FunctionType

     
def text (filename ='data.txt'):
    Material=Alloy()
    R=Factors.R
    Q1=Factors.Q1
    Q2=Factors.Q2
    Q3=Factors.Q3
    Composition=Material.w
    
    os.chdir('../data') # moving to the directory data
    # os.chdir('../scripts') # moving to the directory scripts
    fh = open (filename, 'w')
    fh.write('//********************************************************** \n')
    fh.write('//***Coefficients necessary to define model in FlowVision*** \n')
    fh.write('//********************************************************** \n')
    fh.write('\n')
    fh.write('gs=%d; // Grain size coefficient \n' %Material.gs)
    fh.write('R=%.6f;   // Gas constant\n' %R)
    fh.write('Q1=%.6f;   // Activation energy \n' %Q1)
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
    fh.write('\n')
    fh.write('//Function type - %s \n' % FunctionType)
    fh.write('\n')
    fh.write('//***Chenical composition*** \n')
    fh.write('// ')  
    print(Composition, file=fh)
      
    fh.close()
    return fh

def results2text(data, filename ='results.txt'):
    fh = open (filename, 'w')
    print(data.iloc[-1, :], file=fh)
    fh.write('\n')
    fh.write('//***Chenical composition*** \n')
    fh.write('// ')
    print(Composition, file=fh)
    fh.write('\n')
    fh.write('Austenisation Temperature = %.2f; \n' % T[0])
    fh.write('Cooling time = %.2f; \n' % t[1])
    Cooling_rate=(T[0]-T[1])/t[1]
    fh.write('Cooling rate = %.2f; \n' % Cooling_rate)
    fh.close()
    return fh

def TTT2PandasDataframe():
    Steel=Factors()
    # Temperatures
    # Ferrite start and end
    temp_f_start=Steel.temp_f[0]
    temp_f_end=Steel.temp_f[-1]
    temp_f=Steel.temp_f
    # Pearlite start and end
    temp_p_start=Steel.temp_p[0]
    temp_p_end=Steel.temp_p[-1]
    temp_p=Steel.temp_p
    # Bainite start and end
    temp_b_start=Steel.temp_b[0]
    temp_b_end=Steel.temp_b[-1]
    temp_b=Steel.temp_b
     
    # Start and end temperatures in the TTT diagramm
    temp_start=Steel.temp_f[0]
    temp_end=Steel.temp_b[-1]
    Temperature=np.arange(temp_start, temp_end, -delta_T)
    
    # Start Curves
    start_f=Steel.tau_f[:, 0]
    start_p=Steel.tau_p[:, 0]
    start_b=Steel.tau_b[:, 0]
         
    # End Curves
    end_f=Steel.tau_f[:, -1]
    end_p=Steel.tau_p[:, -1]
    end_b=Steel.tau_b[:, -1]
    
    # find indexes for each phase
    index_f=np.arange(0, len(temp_f), 1)
    index_b=np.arange(index_f[-1]+1, len(Temperature), 1)
    index_p=np.arange(len(temp_f)-len(temp_p), len(temp_f), 1)
    index_TTT=np.concatenate((index_f,index_b), axis=None)
    
       
    # create Pandas dataframe to collect data
    TTT = pd.DataFrame(columns=['Temperature', 
                              'FerriteStart','FerriteEnd', 
                              'PearliteStart','PearliteEnd',
                              'BainiteStart', 'BainiteEnd'], index=index_TTT)
    # write Temperature array   
    TTT['Temperature']=Temperature
    
    # write Ferrite start and end curves
    for i in index_f:
        TTT.loc[i,'FerriteStart']=start_f[i]
        TTT.loc[i,'FerriteEnd']=end_f[i]
        
    # write Pearlite start and end curves
    j=-1    
    for i in index_p:
        j+=1
        TTT.loc[i,'PearliteStart']=start_p[j]
        TTT.loc[i,'PearliteEnd']=end_p[j]
        
    # write Bainite start and end curves    
    j=-1 
    for i in index_b:
        j+=1
        TTT.loc[i,'BainiteStart']=start_b[j]
        TTT.loc[i,'BainiteEnd']=end_b[j]
        
    return TTT

def TTT2Excel():
    TTT=TTT2PandasDataframe()
    os.chdir('../data') # moving to the directory data
    TTT.to_excel("TTT.xlsx")
    

def tau_n_b_2Excel (dataframe, filename, sheet_name):
        
    os.chdir('../data') # moving to the directory data
    
       
    Temp_k_n__b_feq=pd.DataFrame(data=[dataframe.temp, dataframe.tau, 
                                    dataframe.n,dataframe.b, dataframe.f_eq]).transpose()
 


    
    

    Temp_k_n__b_feq.to_excel(filename, sheet_name=sheet_name, index=False)
             
    return Temp_k_n__b_feq
            
def tau_2text(dataframe, filename):
    os.chdir('../data') # moving to the directory data
    
    Temp_tau=pd.DataFrame(data=[dataframe.temp, dataframe.tau]).transpose()
    Temp_tau.to_csv(filename, header=None, index=None, sep=' ', mode='w')
    
       
    return #Temp_tau


def n_2text(dataframe, filename):
    os.chdir('../data') # moving to the directory data
    
    Temp_n=pd.DataFrame(data=[dataframe.temp, dataframe.n]).transpose()
    Temp_n.to_csv(filename, header=None, index=None, sep=' ', mode='w')
    
    
    
    return #Temp_n

def b_2text(dataframe, filename):
    os.chdir('../data') # moving to the directory data
    
    Temp_b=pd.DataFrame(data=[dataframe.temp, dataframe.b]).transpose()
    Temp_b.to_csv(filename, header=None, index=None, sep=' ', mode='w')
       
    return #Temp_b

def time_s_2text(dataframe, filename):
    os.chdir('../data') # moving to the directory data
    
    Temp_time_s=pd.DataFrame(data=[dataframe.temp, dataframe.time_s]).transpose()
    Temp_time_s.to_csv(filename, header=None, index=None, sep=' ', mode='w')
       
    return #Temp_time_s




def eq_2text(dataframe, filename):
    os.chdir('../data') # moving to the directory data
    
    Temp_f_eq=pd.DataFrame(data=[dataframe.temp, dataframe.f_eq]).transpose()
    Temp_f_eq.to_csv(filename, header=None, index=None, sep=' ', mode='w')
    
    return #Temp_f_eq

def phasa_points(data, filename ='Ae3_Ae1_Bs_Ms.txt'):
    os.chdir('../data') # moving to the directory data
    
    # calculating necessary constants
    Ae3=data.Ferrite.temp[0]
    Ae1=data.Pearlite.temp[0]
    Bs=data.Bainite.temp[0]
    Ms=data.Bainite.iloc[-1, 0]

    # writing everything into text file
    fh = open (filename, 'w')
    fh.write('//********************************************************** \n')
    fh.write('//***Coefficients necessary to define model in FlowVision*** \n')
    fh.write('//********************************************************** \n')
    fh.write('\n')
    fh.write('Ae3=%.2f; // Ferrite start temperature \n' %Ae3)
    fh.write('Ae1=%.2f;   // Pearlite start temperature\n' %Ae1)
    fh.write('Bs=%.2f;   // Bainite start temperature \n' %Bs)
    fh.write('Ms=%.2f;   // Martensite start temperature \n' %Ms)
    fh.close()
    return fh








           
  
    
  
    
  
    
  