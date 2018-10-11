from .SymEngineFast import *
from . import Diag_Schar
import symengine as si
from . import PyxSchar
from pyx import *
import numpy as np
import os
'''
I wrote this programm to plot the Action against one variable of choice.
Here it's the impuls, it's either K1 or K2.
Maybe one next step would be to generalize this programm,
so i can plot the action also against the weigts or frequenz.

Also i did this only for N = 1 and 2.

'''
def func(K_List, Rho_List, w_List, kappa,k1):
    K_List[ind]= k1
    print('K_List=', K_List, 'WertListe=',WertListe, 'yaaaaaaay1')
    #Wert = get_Action(N, Integration_bound, T, K_List, Rho_List,w_List,kappa,
    #        False, False, 1)
    CFS_Action.K_Liste = K_List
    CFS_Action.Rho_Liste = Rho_List
    CFS_Action.w_Liste = w_List
    CFS_Action.kappa = kappa 
 
    Wert = CFS_Action.get_Action()
    print(Wert)
    WertListe.append(Wert)
    KK.append(k1)
    return WertListe, KK
    ###return WertListe, KK

if __name__ == "__main__":
    T = 2*np.pi             #float, Lifetime of the universe, which is needed for the
                            #       Schwartzfunktion
    N = 2                   #integer, Shell-Number

    kappa = 0.0001            #float, Needed for boundary constant
    kappa_Anzahl = 1        #float, Needed for plotting a diagram for various
                            # kappas

    x_Anf = 0
    x_End = np.pi
    Integration_bound = [[x_Anf, x_End], [0,2*np.pi]]

    ind = 0                 #Which Impuls should be varied, 0 for K1 and 1 for K2


    w_List = [0.,1.]
    K_List = [0.,0.]
    Number_of_Weights = 11
    Kurve_Names = []
    PDFTitle = 'Rout_Wechs_N_%d_GewichteScharAnz%dVarK_%1dKappa_%1.5f'%(N,Number_of_Weights,
            ind+1, kappa)
    Rho_List = [1,0] # i have to set this here for the initialing of the class C_FS
    c = PyxSchar.initialXY(0,20,0,4, r'$K_%1d$'%(ind+1),'Action',10,10 , 'tr')
    System_Parameters = [K_List, Rho_List, w_List, kappa]
    CFS_Action = C_F_S(N, Integration_bound, T, System_Parameters, Schwartzfunktion = True, 
                 Comp_String = 1, Integration_Type = 1)
    
    dd=[]

    for k in range(Number_of_Weights):
        rho1 = 0.1*k
        Rho_List =[rho1, (1- rho1)/3]
        dk = 0.1
        k1 = 0.
        WertListe = []
        KK =[]
        while 1:
            WertListe, KK =func(K_List, Rho_List, w_List, kappa, k1)
            print('WertListe', WertListe)
            l1 = len(WertListe) - 1
            if  l1>=2 :
                a = abs((WertListe[l1-2] - WertListe[l1-1])/(KK[l1-2] -
                    KK[l1-1]))
                b = abs((WertListe[l1] - WertListe[l1-1])/(KK[l1] - KK[l1-1]))
                if  b != 0 and (a/b < 0.5 or a/b >2):
                    WertListe.pop(l1)
                    KK.pop(l1)
                    k1 -= dk
                    dk = 0.02
                    k1 += dk
                    WertListe, KK = func(K_List, Rho_List, w_List, kappa, k1)
                    print(1,dk)
                elif a != 0 and (b/a < 0.5 or b/a > 2):
                    WertListe.pop(l1)
                    KK.pop(l1)
                    k1 -= dk
                    dk = 0.02
                    k1 += dk
                    WertListe, KK = func(K_List, Rho_List, w_List, kappa, k1)
                    print(2,dk)
                else:
                    dk = dk*2
                    dk = min(dk, 0.5)
                    print(3,dk)
            if WertListe[l1] > 1.5 or k1 > 20:
                break

            k1 += dk
        LL = np.array(WertListe)
        LL2 = LL/LL[0]
        print('LL=', LL,'KK=', KK)
        dd.append(PyxSchar.set_values(KK, LL2,
            r"$\rho_1=%2.4f, \rho_2=%2.4f$"%(Rho_List[0], Rho_List[1])))

    
    if ind ==0:
        PyxSchar.plot_diag(dd, '$kappa=%2.5f, K_{%1d}=%2.2f$'%(kappa, 2, 0),(5,10), PDFTitle, c)
    else:
        PyxSchar.plot_diag(dd, '$kappa=%2.5f, K_{%1d}=%2.2f$'%(kappa, 1, 0),
            (5,10), PDFTitle, c)

    f = open(PDFTitle+'.pdf', 'a')
    f.write('x_values'+str(KK)+'\n')
    f.write('y_values'+ str(LL2)+'\n'+'unnormed values\n' +str(LL))
    f.close()
