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
       #self._Ferrite=self.get_variables()[0]
       #self._Pearlite=self.get_variables()[1]
       #self._Bainite=self.get_variables()[2]
       ######################################
       ######################################
       
       self.Ferrite=self.k_n_Error()[0]
       self.Pearlite=self.k_n_Error()[1]
       self.Bainite=self.k_n_Error()[2]
       #self.results=self.k_n_Error()[3]

       
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
        Temperature = a
        temp_f_s, time_f_s, temp_f_e, time_f_e = b, c, d, e
        temp_p_s, time_p_s, temp_p_e, time_p_e = f, g, h, i
        temp_b_s, time_b_s, temp_b_e, time_b_e = j, k, l, m
        
        # add points for better sollution when we define k and n in 
        # JMAK equation
        # 10% point
        coef1=0.25
        
        time_f_01=coef1*(time_f_e-time_f_s)
        time_p_01=coef1*(time_p_e-time_p_s)
        time_b_01=coef1*(time_b_e-time_b_s)
        
        # 50% point
        time_f_05=0.5*(time_f_e-time_f_s)
        time_p_05=0.5*(time_p_e-time_p_s)
        time_b_05=0.5*(time_b_e-time_b_s)
        
        # 90% point
        coef2=0.75
        time_f_09=coef2*(time_f_e-time_f_s)
        time_p_09=coef2*(time_p_e-time_p_s)
        time_b_09=coef2*(time_b_e-time_b_s)    
        
     
         
        ######################################################################
        # Create dataframe with values for Ferrite 
        Ferrite=pd.DataFrame(data=[temp_f_s, time_f_s, time_f_01, time_f_05,
                                   time_f_09, time_f_e]).transpose()
        Ferrite.columns = ['temp', 'time_s', 'time_01', 'time_05', 'time_09','time_e']
        
        ###################
        #creating columns for k, n and Error for the further calculating of
        # these coefficients
        Ferrite['k']=np.zeros(len(Ferrite['temp']))
        Ferrite['n']=np.zeros(len(Ferrite['temp']))
        Ferrite['Error']=np.zeros(len(Ferrite['temp']))
        Ferrite['f_eq'] = np.ones(len(Ferrite['temp']), dtype=float)

        ######################################################################
        # Create dataframe with values for Pearlite
        Pearlite=pd.DataFrame(data=[temp_p_s, time_p_s, time_p_01, time_p_05, 
                                    time_p_09,time_p_e]).transpose()
        Pearlite.columns = ['temp', 'time_s', 'time_01', 'time_05', 'time_09','time_e']
        #####
        Pearlite['k']=np.zeros(len(Pearlite['temp']))
        Pearlite['n']=np.zeros(len(Pearlite['temp']))
        Pearlite['Error']=np.zeros(len(Pearlite['temp']))
        Pearlite['f_eq'] = np.ones(len(Pearlite['temp']), dtype=float)
        
        ######################################################################
        # Create dataframe with values for Bainite
        Bainite=pd.DataFrame(data=[temp_b_s, time_b_s, time_b_01, time_b_05, 
                                   time_b_09, time_b_e]).transpose()
        Bainite.columns = ['temp', 'time_s', 'time_01', 'time_05', 'time_09','time_e']    
        #####
        Bainite['k']=np.zeros(len(Bainite['temp']))
        Bainite['n']=np.zeros(len(Bainite['temp']))
        Bainite['Error']=np.zeros(len(Bainite['temp']))
        Bainite['f_eq'] = np.ones(len(Bainite['temp']), dtype=float)   
                
        return Ferrite, Pearlite, Bainite
    
    
    def k_n_Error(self):
        results=[]
        for i in range (3):
            dataframe = self.get_variables()[i]
            t001=dataframe.time_s.to_numpy()
            t099=dataframe.time_e.to_numpy()
            j=len(t001)-1
            dataframe['k'].loc[j] = JMAK.Fit_k_n(t001[j], t099[j])[0]
            dataframe['n'].loc[j] = JMAK.Fit_k_n(t001[j], t099[j])[1]
            dataframe['Error'].loc[j] = JMAK.Fit_k_n(t001[j], t099[j])[2]
            while j >=1:
                j=j-1
                jj=0
                guess=[dataframe['k'].loc[j+1], dataframe['n'].loc[j+1]]
                k, n, Error = JMAK.Fit_k_n(t001[j], t099[j], guess)[0:3]
                ## If we can't find the solution we increase the k guess value
                ## for 1000**jj (1000**1, 1000**2, ..., 1000**5)
                if Error<=0 or Error > 0.5:
                    while jj<6:
                        jj=jj+1
                        guess=[dataframe['k'].loc[j+1]*1000**jj, dataframe['n'].loc[j+1]]
                        k, n, Error = JMAK.Fit_k_n(t001[j], t099[j], guess)[0:3]
                        if Error>0 and Error < 0.5: jj=6
                                           
                dataframe['k'].loc[j]=k
                dataframe['n'].loc[j]=n
                dataframe['Error'].loc[j]=Error
                
            results.append(dataframe)
        Ferrite=results[0]
        Pearlite=results[1]
        Bainite=results[2]
        return Ferrite, Pearlite, Bainite, results
            

# data=Create_variables_from_Excel ()
# Pearlite=data.Pearlite
# Ferrite=data.Ferrite
# Bainite=data.Bainite

#Pearlite = Create_variables_from_Excel().Pearlite





