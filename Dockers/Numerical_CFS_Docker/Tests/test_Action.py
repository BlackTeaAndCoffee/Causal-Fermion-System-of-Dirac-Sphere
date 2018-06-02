from Numerical_CFS_Docker.configfunktion import writeconfig
from Numerical_CFS_Docker.SymEngineFast import *
import numpy as np

def func(K_List, Rho_List, w_List, kappa, k1, ind, CFS_Action):
    K_List[ind]= k1
    
    CFS_Action.K_Liste = K_List
    CFS_Action.Rho_Liste = Rho_List
    CFS_Action.w_Liste = w_List
    CFS_Action.kappa = kappa 
 
    Wert = CFS_Action.get_Action()
    return Wert 

def main():
    T = 1                   #float, Lifetime of the universe, which is needed for the
                            #       Schwartzfunktion
    N = 2                   #integer, Shell-Number

    kappa = 0            #float, Needed for boundary constant
    kappa_Anzahl = 1        #float, Needed for plotting a diagram for various
                            # kappas

    x_Anf = 0
    x_End = np.pi
    Integration_bound = [[x_Anf, x_End], [0,T]]

    ind = 0                 #Which Impuls should be varied, 0 for K1 and 1 for K2


    w_List = [0.,1.]
    K_List = [0.,0.]
    Rho_List = [1,0] # i have to set this here for the initialing of the class C_FS
    System_Parameters = [K_List, Rho_List, w_List, kappa]
    CFS_Action = C_F_S(N, T, System_Parameters,Integration_bound,  Schwartzfunktion = False, 
                 Comp_String = False, Integration_Type = 1)
    
    dd=[]

    WertListe = []
    for k in range(20):
        k1 = k/10
        KK =[]
        WertListe.append(func(K_List, Rho_List, w_List, kappa, k1, ind, CFS_Action) + 1 - 1/(2*np.pi))
    print('WertListe', WertListe)
    return WertListe

def test_action():
    assert main() == [0.0008277342004438248, 0.0008309955105414712, 0.000840194562103519, 
                      0.0008536589176576792, 0.0008689658481976503, 0.0008835306967854006, 
                      0.0008951221231393858, 0.0009022625399987272, 0.0009042350910075936, 
                      0.0009010571906092057, 0.0008931757689961961, 0.0008812958950642769, 
                      0.0008662152758828479, 0.000848700978669209, 0.000829473187787505, 
                      0.0008091024283588821, 0.0007880785363976461, 0.0007667879751987183, 
                      0.0007455241133152157, 0.0007245150453482341]

