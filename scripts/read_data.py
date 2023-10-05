#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 19:22:49 2023

@author: aas
"""


import os.path
import pandas as pd
from solve_JMAK import JMAK
import numpy as np
# from plot import Plot



class Create_variables_from_Excel ():
    
    def __init__(self, filename='TTT.xlsx', filepath='../data'):
       self._filename = filename
       self._filepath = filepath
       ######################################
       ######################################
       self._dataTTT=self.get_TTT()[0]
       self._curve_names=self.get_TTT()[1]
       ######################################
       ######################################
       # self._Ferrite=self.get_variables()[0]
       # self._Pearlite=self.get_variables()[1]
       # self._Bainite=self.get_variables()[2]
       ######################################
       ######################################
       
       self.Ferrite=self.tau_n_b()[0]
       self.Pearlite=self.tau_n_b()[1]
       self.Bainite=self.tau_n_b()[2]
       self.results=self.tau_n_b()[3]

       
    def get_TTT(self):
       """
       create dataframe from Excel file and a list of curve names
       """
       os.chdir(self._filepath) # moving to the directory data
       data_TTT=pd.read_excel(self._filename, index_col=0) #create dataframe with data
       curve_names = data_TTT.columns.values.tolist()
       return data_TTT, curve_names
   
    def get_variables(self):
        """
        create a separate dataframe for each phase from the common one
        """
        data_TTT=self._dataTTT
        curve_names=self._curve_names
        data={}
        data[curve_names[0]]=data_TTT[curve_names[0]].to_numpy(dtype=float)
        for i in range (1, len(curve_names)):
            name=curve_names[i]
            dataframe=pd.DataFrame (data = [data[curve_names[0]], 
                                            data_TTT[name]]).transpose()
            dataframe=dataframe.dropna(axis=0)
            temp=dataframe.to_numpy(dtype=float)[:, 0]
            time=dataframe.to_numpy(dtype=float)[:, 1]
            data[name + ' temp']=temp
            data[name + ' time']=time
            
        # split the dictionary to separate variables
        # by default we assume that we have Ferrite, Pearlite, Bainite
        # the "resе" variable for the case we have something in addition
        a, b, c, d, e, f, g, h, i, j, k, l, m, *rest = data.values()
        # Temperature = a
        temp_f_s, time_f_s, temp_f_e, time_f_e = b, c, d, e
        temp_p_s, time_p_s, temp_p_e, time_p_e = f, g, h, i
        temp_b_s, time_b_s, temp_b_e, time_b_e = j, k, l, m
        
        # middle point
        time_f_05=0.5*(time_f_e-time_f_s)
        time_p_05=0.5*(time_p_e-time_p_s)
        time_b_05=0.5*(time_b_e-time_b_s)
        
                
        ######################################################################
        # Create dataframe with values for Ferrite 
        Ferrite=pd.DataFrame(data=[temp_f_s, time_f_s, time_f_05, time_f_e]).transpose()
        Ferrite.columns = ['temp', 'time_s', 'time_05','time_e']
        ##################
        Ferrite['X001']=np.zeros(len(Ferrite['temp']))
        Ferrite['X050']=np.zeros(len(Ferrite['temp']))
        Ferrite['X099']=np.zeros(len(Ferrite['temp']))
               
        ###################
        #creating columns for tau, n and b for the further calculating of
        # these coefficients
        Ferrite['tau']=np.zeros(len(Ferrite['temp']))
        Ferrite['n']=np.zeros(len(Ferrite['temp']))
        Ferrite['b']=np.zeros(len(Ferrite['temp']))
        Ferrite['f_eq'] = np.ones(len(Ferrite['temp']), dtype=float)

        ######################################################################
        # Create dataframe with values for Pearlite
        Pearlite=pd.DataFrame(data=[temp_p_s, time_p_s, time_p_05,time_p_e]).transpose()
        Pearlite.columns = ['temp', 'time_s', 'time_05', 'time_e']
        #####
        Pearlite['X001']=np.zeros(len(Pearlite['temp']))
        Pearlite['X050']=np.zeros(len(Pearlite['temp']))
        Pearlite['X099']=np.zeros(len(Pearlite['temp']))
        ######
        Pearlite['tau']=np.zeros(len(Pearlite['temp']))
        Pearlite['n']=np.zeros(len(Pearlite['temp']))
        Pearlite['b']=np.zeros(len(Pearlite['temp']))
        Pearlite['f_eq'] = np.ones(len(Pearlite['temp']), dtype=float)
        
        ######################################################################
        # Create dataframe with values for Bainite
        Bainite=pd.DataFrame(data=[temp_b_s, time_b_s, time_b_05, time_b_e]).transpose()
        Bainite.columns = ['temp', 'time_s','time_05','time_e']    
        #####
        Bainite['X001']=np.zeros(len(Bainite['temp']))
        Bainite['X050']=np.zeros(len(Bainite['temp']))
        Bainite['X099']=np.zeros(len(Bainite['temp']))
                #####
        Bainite['tau']=np.zeros(len(Bainite['temp']))
        Bainite['n']=np.zeros(len(Bainite['temp']))
        Bainite['b']=np.zeros(len(Bainite['temp']))
        Bainite['f_eq'] = np.ones(len(Bainite['temp']), dtype=float)   
                
        return Ferrite, Pearlite, Bainite
    
    


    def tau_n_b(self):
        results=[]
        for i in range (3):
            dataframe = self.get_variables()[i]
            t001=dataframe.time_s.to_numpy()
            t099=dataframe.time_e.to_numpy()
            t050=dataframe.time_05.to_numpy()
            j=len(t001)
          
            while j >=1:
                j=j-1
                b, n, tau = JMAK.b_n_tau_incub (t001[j], t099[j])
                # X001=JMAK.equation_incubation_tau(t001[j], tau, n, t001[j])
                # X050=JMAK.equation_incubation_tau(t050[j], tau, n, t001[j])
                # X099=JMAK.equation_incubation_tau(t099[j], tau, n, t001[j])
                X001=JMAK.equation_incubation_tau(t001[j], tau, n, t001[j])
                X099=JMAK.equation_incubation_tau(t099[j], tau, n, t001[j])
                X050=JMAK.equation_incubation_tau(t050[j], tau, n, t001[j])
                
                
                
                
                dataframe['tau'].loc[j]=tau
                dataframe['n'].loc[j]=n
                dataframe['b'].loc[j]=b
                dataframe['X001'].loc[j]=X001
                dataframe['X050'].loc[j]=X050
                dataframe['X099'].loc[j]=X099
              
            results.append(dataframe)
        Ferrite=results[0]
        Pearlite=results[1]
        Bainite=results[2]
        return Ferrite, Pearlite, Bainite, results








# #%%

# data=Create_variables_from_Excel ()
# Pearlite=data.Pearlite
# Ferrite=data.Ferrite
# Bainite=data.Bainite

# #%%
# i=0
# t001=Ferrite.time_s[i]
# t099=Ferrite.time_e[i]
# ratio=t099/t001
# Temp=Ferrite.temp[i]

# b, n, tau = JMAK.b_n_tau_incub(t001, t099)

# #time=JMAK.x_t_interpalation(t001, t099)[0]
# time=np.linspace(t001, t099)
# X=JMAK.equation_incubation_tau(time, tau, n, t001)
# figure=Plot.plot_isothermal_transformation(time, X, Temp)


