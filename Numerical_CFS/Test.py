from sympy.printing import ccode
from scipy import integrate
from symengine import I
from sympy import *
from Numerical_CFS.configfunktion import configfunktion
from Numerical_CFS.Lib_Action_Minimum import *
import configparser
import sympy as sy
import symengine as si
import numpy as np
import random
import time
import os
import sys
import scipy.misc
import ctypes
from Numerical_CFS.SymEngineFast import *

if __name__ == "__main__":
    var_K, var_Rho, var_w = configfunktion('Vary_Parameters') #boolean
    K_Anf, K_End, pre_K_List= configfunktion('Impuls') # floats and List
    w_Anf, w_End, pre_w_List= configfunktion('Frequenz') # flaots and List
    Constant, kappa, pre_Rho_List = configfunktion('Constraints') #floats and List
    Anzahl_N, first = configfunktion('System_sizes') # integer and integer
    StartWithGivenMinima, Test_Action = configfunktion('Test') #boolean, boolean
    random_K, random_Rho, random_w =  configfunktion('Set_Initialstate_randomly')# boolean


   
    delta_K = 5/10
    delta_w = 1/10
    delta_Rho = 1/10
    
    T = 2*np.pi

    x_Anf = 0
    x_End = np.pi

    Integration_bound = [[x_Anf, x_End], [0,2*np.pi]]
    Wirk = []
    Boltzmann_Constant =0.00005
    Mittelgr = 4
    for_rho = 1

    SN = 5
    #Minimum =[[[0.61402128722400451, 5.4919577589099093], [0.96197069,  0.01267644], [0, 1]], 0.007515003816659801]
    Constant = 1
#    pre2_w_List = eval(pre_w_List)

    K_List = eval(pre_K_List)
    K_List[4]=1
    Rho_Values = Rho_Class(SN, Constant)
    w_List = eval(pre_w_List) # I put this list assignment here,
                                    #I set the list in settings.cfs, and it's like 
                                    #[i for i in range(SN)]. So i need SN.
    #Rho_List = Rho_Values(np.array([0.31848333, 0.21734017, 0.00491603, 0.00001, 0.000002, 0.0000001]))
    Rho_List = Rho_Values(np.array([0.29522226, 0.22771744, 0.00360424, 0.0001, 0.000000005]))
    #Rho_List = Rho_Values.SameRhoValues()
    System_Parameters = [K_List, Rho_List, w_List, kappa]
    print(System_Parameters)

    CFS_Action = C_F_S(SN, Integration_bound, T, System_Parameters, Schwartzfunktion = True, 
    Comp_String = False, Integration_Type = 1, Test_Action = False)
    CFS_Action.get_Action()
  



