#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 15:32:22 2023

@author: aas
"""
import numpy as np
from sigmoidal import Sigmoidal
from initial_data import Composition, t, T, delta_x, delta_T, delta_t

#%%

class Alloy(Sigmoidal):
    """
    Alloy properties (composition in wt.% and prior austenite grain size)
    """
    def __init__(self, w=Composition):
        
        self.Sigmoidal=Sigmoidal()
        # Alloy composition and n coefficients
        self.w = w
        self.n1_P=0.32
        self.n1_F=0.41
        self.n1_B=0.29
        self.n2_P=3
        self.n2_F=3
        self.n2_B=2
        # Main elements
        self.C = w.get('C', 0)
        self.Mn = w.get('Mn', 0)
        self.Si = w.get('Si', 0)
        self.Ni = w.get('Ni', 0)
        self.Cr = w.get('Cr', 0)
        self.Mo = w.get('Mo', 0)
        self.Co = w.get('Co', 0)
        self.As = w.get('As', 0)
        self.W = w.get('As', 0)
        self.V = w.get('V', 0)
        self.Cu=w.get('Cu', 0)
        # Grain size
        self.gs = w.get('GrainSize', 0)
        #Calculated coeficients
        
        # Check if we use S(x) or I(x) function (see module Sigmoidal for details)
        if self.Sigmoidal.f==self.Sigmoidal.S:
            self.FC=self.FC_S()
            self.PC=self.PC_S()
            self.BC=self.BC_S()
        if  self.Sigmoidal.f==self.Sigmoidal.I:
            self.FC=self.FC_I()
            self.PC=self.PC_I()
            self.BC=self.BC_I()
         
        self.alpha_martensite = self.alpha_martensite_VanBohemen()
        #Temperatures Start/End
        self.Ae3=self.Ae3_Andrews()
        self.Ae1=self.Ae1_Andrews()
        self.Bs=self.Bs_Li()
        self.Ms=self.Ms_VanBohemen()
        
    @staticmethod 
    def GrainSize(w):
        return w.get('GrainSize', 0)
        
        
    def FahrenheitToCelsius(self, TF):
        """
        Converts temperature in Fahrenheit to Celsius
        """
        return (TF - 32.)*5./9.


    def FC_S(self):
        """
        Function that expresses the effects of the alloying elements on
        on the kinetics of ferrite transformation (Li aproach)
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return np.exp((1.0 + 6.31*C + 1.78*Mn + 0.31*Si + 
                       1.12*Ni + 2.7*Cr + 4.06*Mo))

    def FC_I(self):
        """
        Function that expresses the effects of the alloying elements on
        on the kinetics of ferrite transformation (Saunders aproach)
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return 60.*Mn+2.*Ni+68.*Cr+244.*Mo



    def PC_S(self):
        """
        Function that expresses the effects of the alloying elements on
        on the kinetics of pearlite transformation
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return np.exp(-4.25 + 4.12*C + 4.36*Mn + 0.44*Si + 1.71*Ni + 3.33*Cr + 5.19*np.sqrt(Mo))

    def PC_I(self):
        """
        Function that expresses the effects of the alloying elements on
        on the kinetics of perlite transformation (Saunders aproach)
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return 1.8+5.4*(Cr+Mo+4.*Mo*Ni)

    

    def BC_S(self):
        """
        Function that expresses the effects of the alloying elements on
        on the kinetics of bainite transformation
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return np.exp(-10.23 + 10.18*C + 0.85*Mn + 0.55*Ni + 0.9*Cr + 0.36*Mo)

    def BC_I(self):
        """
        Function that expresses the effects of the alloying elements on
        on the kinetics of bainite transformation (Saunders aproach)
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return (2.3 + 10.*C + 4.*Cr + 19.*Mo)*10**(-4)


    def Ae1_Grange(self):
        """
        Grange's equation for Ae1
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return self.FahrenheitToCelsius(1333 - 25*Mn + 40*Si - 26*Ni + 42*Cr)


    def Ae3_Grange(self):
        """
        Grange's equation for Ae3
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        self.FahrenheitToCelsius(1570 - 323*C - 25*Mn + 80*Si - 32*Ni - 3*Cr)
        return self.FahrenheitToCelsius(1570 - 323*C - 25*Mn + 80*Si - 32*Ni - 3*Cr)


    def Ae1_Andrews(self):
        """
        Andrews' equation for Ae1
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return 723 - 16.9*Ni + 29.1*Si + 6.38*W - 10.7*Mn + 16.9*Cr + 290*As


    def Ae3_Andrews(self):
        """
        Andrews' equation for Ae3
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        V=self.V
        Cu=self.Cu
        gs=self.gs
        return 910 - 203*np.sqrt(C) + 44.7*Si - 15.2*Ni + 31.5*Mo + 104*V + 13.1*W - \
            30.0*Mn + 11.0*Cr + 20.0*Cu  # - 700*P - 400*Al - 120*As - 400*Ti


    def Bs_Li(self):
        """
        Bainite start calculation from Li's work
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return 637 - 58*C - 35*Mn - 15*Ni - 34*Cr - 41*Mo
    

    def Bs_VanBohemen(self):
        """
        Bainite start calculation from Van Bohemen's work
        [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return 839 - (86*Mn + 23*Si + 67*Cr + 33*Ni + 75*Mo) - 270*(1 - np.exp(-1.33*C))


    def Ms_Andrews(self):
        """
        Andrews' equation for Ms
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        
        return 539 - 423*C - 30.4*Mn - 17.7*Ni - 12.1*Cr - 7.5*Mo + 10*Co - 7.5*Si


    def alpha_martensite_VanBohemen(self):
        """
        Martensite transformation rate constant
        [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        
        return 1e-3*(27.2 - (0.14*Mn + 0.21*Si + 0.11*Cr + 0.08*Ni + 0.05*Mo) - 19.8*(1-np.exp(-1.56*C)))


    def Ms_VanBohemen(self):
        """
        Martensite start temperature
        [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
        """
        C=self.C
        Mn=self.Mn
        Si=self.Si
        Ni=self.Ni
        Cr=self.Cr
        Mo=self.Mo
        Co=self.Co
        As=self.As
        W=self.W
        gs=self.gs
        return 565 - (31*Mn + 13*Si + 10*Cr + 18*Ni + 12*Mo) - 600*(1-np.exp(-0.96*C))
    

    # def Hv_martensite(phi700, **comp):
    #     """
    #     Martensite Vickers hardness empirical equation
    #     (Maynier et al.)
    #     """
    #     C = comp.get('C', 0)
    #     Mn = comp.get('Mn', 0)
    #     Si = comp.get('Si', 0)
    #     Ni = comp.get('Ni', 0)
    #     Cr = comp.get('Cr', 0)
    #     return 127 + 949*C + 27*Si + 11*Mn + 8*Ni + 16*Cr + 21*np.log10(phi700*3600)


    # def Hv_bainite(phi700, **comp):
    #     """
    #     Bainite Vickers hardness empirical equation
    #     (Maynier et al.)
    #     """
    #     C = comp.get('C', 0)
    #     Mn = comp.get('Mn', 0)
    #     Si = comp.get('Si', 0)
    #     Ni = comp.get('Ni', 0)
    #     Cr = comp.get('Cr', 0)
    #     Mo = comp.get('Mo', 0)
    #     return -323 + 185*C + 330*Si + 153*Mn + 65*Ni + 144*Cr + 191*Mo + \
    #         (89 + 53*C - 55*Si - 22*Mn - 10*Ni - 20*Cr - 33*Mo)*np.log10(phi700*3600)


    # def Hv_ferrite_pearlite(phi700, **comp):
    #     """
    #     Ferrite + pearlite Vickers hardness empirical equation
    #     (Maynier et al.)
    #     """
    #     C = comp.get('C', 0)
    #     Mn = comp.get('Mn', 0)
    #     Si = comp.get('Si', 0)
    #     Ni = comp.get('Ni', 0)
    #     Cr = comp.get('Cr', 0)
    #     Mo = comp.get('Mo', 0)
    #     V = comp.get('V', 0)
    #     return 42 + 223*C + 53*Si + 30*Mn + 12.6*Ni + 7*Cr + 19*Mo + \
    #         (10 - 19*Si + 4*Ni + 8*Cr + 130*V)*np.log10(phi700*3600)
    
