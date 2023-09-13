#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 15:34:55 2023

@author: aas
"""


import numpy as np
import pandas as pd
from read_data import Create_variables_from_Excel as create
from solve_JMAK import JMAK


class Equilibria (create):
    
    '''
    Functions that caclulates equilibria fractions of Ferrite and Pearlite at
    temperatures where both phases can rise simultaniously
    
    '''
    
  
    def __init__(self, filename='TTT.xlsx', filepath='../data'):
        
        
        self.data=create(filename, filepath)
        # self.Ferrite = self.data.Ferrite
        # self.Pearlite = self.data.Pearlite
       
        
        ##################################
        # self.temp_intersection = self.get_intersaction()[0]
        # self.f_ = self.get_intersaction()[0]
        # self.p_ = self.get_intersaction()[1]
        ##################################
        ##################################
        self.f_=self.get_equilibria()[0]
        self.p_=self.get_equilibria()[1]
        
        #################################
        #################################
        self.Ferrite=self.fill_eq()[0]
        self.Pearlite=self.fill_eq()[1]
        self.Bainite=self.data.Bainite
        
        

    def get_intersaction(self):
        
        '''
        Search for temperatures where Ferrite and Pearlite transformations can
        happen simultaniously
        
        Input:
        Ferrite and Pearlite dataframes
        
        Output: 
        two dataframes of the same size with Pearlite and Ferrite
        data, including temperaure, time of start and end transformation,
        k and n coefficients
        '''
        
        i_length=max(len(self.data.Ferrite.index), len(self.data.Pearlite.index))
        j_length=min(len(self.data.Ferrite.index), len(self.data.Pearlite.index))
        temp_intersection = np.zeros(j_length)
        f_=np.zeros([j_length, 4])
        p_=np.zeros([j_length, 4])
        for i in range (i_length):
            for j in range (j_length):
                if self.data.Ferrite.temp.loc[i]==self.data.Pearlite.temp.loc[j]:
                    temp_intersection[j]=self.data.Pearlite.temp.loc[j]
                    #############
                    #############
                    f_[j]=[self.data.Ferrite.time_s.loc[i], self.data.Ferrite.time_e.loc[i],
                               self.data.Ferrite.k.loc[i], self.data.Ferrite.n.loc[i]]
                                        
                    p_[j]=[self.data.Pearlite.time_s.loc[j], self.data.Pearlite.time_e.loc[j],
                               self.data.Pearlite.k.loc[j], self.data.Pearlite.n.loc[j]]
        ### put everything into pandas dataframe
        # Ferrite
        f_int=pd.DataFrame(data=f_)
        f_int.columns = ['time_s', 'time_e', 'k', 'n']
        f_int.insert(0, 'temp', temp_intersection)
        
        
        # Pearlite
        p_int=pd.DataFrame(data=p_)
        p_int.columns = ['time_s', 'time_e', 'k', 'n']
        p_int.insert(0, 'temp', temp_intersection)
        
        return f_int, p_int #, temp_intersection, f_, p_
    
    def get_equilibria(self):
        '''
        calculates maximum possible Pearlite and Ferrite fraction at constant
        temperature
        
        Input:
        f_ = self.get_intersaction()[0]
        p_ = self.get_intersaction()[1]
        
        Output:
        two dataframes with Pearlite and Ferrite equilibria fractions in the last
        column of each dataframe
        
        '''
        
        
        f_ = self.get_intersaction()[0]
        p_ = self.get_intersaction()[1]
        
        temp=f_['temp']
        tau001_f, tau099_f = f_['time_s'], f_['time_e']
        tau001_p, tau099_p = p_['time_s'], p_['time_e']
        k1, n1 = f_['k'], f_['n']
        k2, n2 = p_['k'], p_['n']
        f_eq=np.ones(len(temp))
        p_eq=np.ones(len(temp))
        i=-1
        
        while i < len(temp)-1:
            i=i+1
            t_s=min(tau001_f[i], tau001_p[i])
            t_e=min(tau099_f[i], tau099_p[i])
            t_av=t_e-(t_e-t_s)*0.5
            t_eq=t_av
            f_eq[i]=JMAK.equation(t_eq, k1[i], n1[i])
            p_eq[i]=JMAK.equation(t_eq, k2[i], n2[i])
            
            #print(i)
            function = lambda t: JMAK.equation(t, k1[i], n1[i]) + JMAK.equation(t, k2[i], n2[i])
            t_s=min(tau001_f[i], tau001_p[i])
            t_e=min(tau099_f[i], tau099_p[i])
            t_av=t_e-(t_e-t_s)*0.5
            f=function(t_av)
            j=0
            while f <= 0.99 or f > 1:
                if f<= 0.99:
                    j=j+1
                    t_s=t_av
                    t_e=t_e
                    t_av=t_e-(t_e-t_s)*0.5
                    f=function(t_av)
                    t_eq=t_av
                if f>=1:
                    j=j+1
                    t_s=t_s
                    t_e=t_av
                    t_av=t_e-(t_e-t_s)*0.5
                    f=function(t_av)
                    t_eq=t_av
                #print(t_eq)
                f_eq[i]=JMAK.equation(t_eq, k1[i], n1[i])
                f_['f_eq']=f_eq
            #print(f_eq[i])
                p_eq[i]=JMAK.equation(t_eq, k2[i], n2[i])
                p_['f_eq']=p_eq
        return f_, p_
    
    
    def fill_eq (self):
        '''
        fill Pearlite and Ferrite equilibria fractions to the main dataframes 
        for Pearlite and Ferrite

        '''
        #i=Phase[Phase['temp']==Temp].index.values[0]
        
        
        Ferrite=self.data.Ferrite
        Pearlite=self.data.Pearlite
              
        for i in range (len(Ferrite.index)):
            for j in range (len(self.f_.index)):
                if Ferrite.temp[i]==self.f_.temp[j]:
                    Ferrite.f_eq[i]=self.f_.f_eq[j]
        for i in range (len(Pearlite.index)):
            for j in range (len(self.p_.index)):
                if Pearlite.temp[i]==self.p_.temp[j]:
                    Pearlite.f_eq[i]=self.p_.f_eq[j]
        return Ferrite, Pearlite
    
                        

# data=Equilibria()
# Ferrite=data.Ferrite
# Pearlite=data.Pearlite
# Bainite=data.Bainite
# f_=data.f_
# p_=data.p_

