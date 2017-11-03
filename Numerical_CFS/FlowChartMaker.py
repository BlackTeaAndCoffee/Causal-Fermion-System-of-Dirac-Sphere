from pycallgraph import PyCallGraph
from pycallgraph import Config
from pycallgraph.output import GraphvizOutput
from pycallgraph import GlobbingFilter
from .Wirkungs_Minimum import *
import symengine as si
import os

if( __name__ == "__main__"):
    var_K, var_Rho, var_w = configfunktion('Vary_Parameters') #boolean
    K_Anf, K_End, pre_K_List= configfunktion('Impuls') # floats and List
    w_Anf, w_End, pre_w_List= configfunktion('Frequenz') # flaots and List
    Constant, kappa, pre_Rho_List = configfunktion('Constraints') #floats and List
    Anzahl_N, first = configfunktion('System_sizes') # integer and integer
    StartWithGivenMinima, Which_Fitness = configfunktion('Test') #boolean, integer
    random_K, random_Rho, random_w =  configfunktion('Set_Initialstate_randomly')# boolean

    variant = which_variant(random_K, random_Rho, random_w)

    print(StartWithGivenMinima, "variant = ", variant)

    T = 2*np.pi
    r = si.symarray('r', 1)
    t = si.symarray('t', 1)

    x_Anf = 0
    x_End = np.pi

    Integration_bound = [[x_Anf, x_End], [0,2*np.pi]]
    Wirk = []
    Selbstval_K_Boltz =0.00005
    Mittelgr = 4
    for_rho = 1


    #Minimum =[[[0.61402128722400451, 5.4919577589099093], [0.96197069,  0.01267644], [0,     1]], 0.007515003816659801]
    Constant = 1
    #    pre2_w_List = eval(pre_w_List)

    pre2_K_List = eval(pre_K_List)

    pre2_Rho_List = eval(pre_Rho_List)

    config = Config()
    config.trace_filter = GlobbingFilter(include=[
            'Wirkungs_Minimum*'
                ])
    graphviz = GraphvizOutput(output_file='PyCallGraphMinimizer.png')

    with PyCallGraph(output=graphviz, config=config):
        MainProg()
