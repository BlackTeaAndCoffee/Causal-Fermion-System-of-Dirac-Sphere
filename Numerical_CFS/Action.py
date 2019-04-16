from sympy.printing import ccode
from scipy import integrate
from sympy import *
from Numerical_CFS.configfunktion import configfunktion
from Numerical_CFS.Lib_Action_Minimum import *
import configparser
import sympy as sy
import numpy as np
import random
import time
import os
import sys
import scipy.misc
import ctypes
from SymEngineFast import *
#from .LibForSimulAnnealing import *


if __name__ == "__main__":
    var_K, var_Rho, var_w = configfunktion('Vary_Parameters') #boolean
    K_Anf, K_End, pre_K_List= configfunktion('Impuls') # floats and List
    w_Anf, w_End, pre_w_List= configfunktion('Frequenz') # flaots and List
    Constant, kappa, pre_Rho_List = configfunktion('Constraints') #floats and List
    Anzahl_N, first, LifeTime = configfunktion('System_sizes') # integer and integer
    StartWithGivenMinima, Test_Action = configfunktion('Test') #boolean, boolean
    random_K, random_Rho, random_w =  configfunktion('Set_Initialstate_randomly')# boolean



    delta_K = 5/10
    delta_w = 1/10
    delta_Rho = 1/10

    T = np.pi

    x_Anf = 0
    x_End = np.pi

    Integration_bound = [[x_Anf, x_End], [0,2*np.pi]]
    Wirk = []
    Boltzmann_Constant =0.00005
    Mittelgr = 4
    for_rho = 1


    #Minimum =[[[0.61402128722400451, 5.4919577589099093], [0.96197069,  0.01267644], [0, 1]], 0.007515003816659801]
    Constant = 1
#    pre2_w_List = eval(pre_w_List)

    pre2_K_List = eval(pre_K_List)

    pre2_Rho_List = eval(pre_Rho_List)



    Rho_Values = Rho_Class(Anzahl_N, Constant)
    vary = Variation_of_Parameters(var_K, var_Rho, var_w, delta_K, delta_Rho, delta_w, Rho_Values)
    pre2_w_List = eval(pre_w_List) # I put this list assignment here,
                                   #I set the list in settings.cfs, and it's like
                                   #[i for i in range(SN)]. So i need SN.
    print('w_List', pre2_w_List)
    Sys_Params= Initial_System_Params(random_K, random_Rho, random_w,
            pre2_K_List, pre2_Rho_List, pre2_w_List, kappa, Constant)

    variant = Sys_Params.variant

    print(StartWithGivenMinima, "variant = ", variant)



    System_Parameters= Sys_Params.Initial_Params_Constructor()

    print('System_Parameters =', System_Parameters)
    CFS_Action = C_F_S(Anzahl_N, T, System_Parameters, Integration_bound, Schwartzfunktion = True,
    Comp_String = 1, Integration_Type = 1, Test_Action = False)


    fitn_wert_y11 = CFS_Action.get_Action()


    System_Parameters[4,0] = fitn_wert_y11
    print('Initial_State', System_Parameters)



