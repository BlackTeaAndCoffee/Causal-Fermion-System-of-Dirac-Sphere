from sympy.printing import ccode
from scipy import integrate
from symengine import I
from sympy import *
from .configfunktion import configfunktion
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
from .SymEngineFast import *
from .LibForSimulAnnealing import *

'''
This Programm contains the code for minimizing the Action.
With the help of some Minimum finder, i want to find ideally
the global Minimum of the Action. The Action as encoded in
SymEngineFast.py has the parameters K_List, Rho_List and w_List.
By parameters i mean, the elements of these Lists.

The Paramater space is in general infinite dimensional. Our idea is
to take N=1 and the corresponding sub-parameter-space and find the minimum
there. Then we expand to the subspace corresponding to N=2, which contains
the subspace for N=1 and find the Minimum there. Ideally we wish
to find some convergency for higher getting N.

The Minimizer we choose is the simulated annealing procedure. The reason
for that is, that the Action must not be smooth.
    * A Problem with that, which already seems to get an issue is that
    each evaluation of the action takes to much time. On Average higher
    N means longer evaluation time. For higher N we also need to probe the
    Paramater space more often. So the simulated annealing method might
    not be good enough.

    *There is another method called bayesian optimization. I haven't tried
    this one out yet. But i will definitely do so soon.

In the simulated annealing algorithm we basically begin with some initial
state. A state is just some combinations of the various paramters in
K_List, Rho_List and w_List.
    The initial state should either be given by the user or the programm should
    pick some state.
        The user would give a state, if he or she want's to give the algorithm
        a starting point. This we will do by each increase of the Shell-Number
        N.

    A state always gets an Energyvalue (float number) and the state with the
    lowest Energy is the desired one. Now we can imagine, that the energy
    landscape has a lot of mountains and valleys and sometime you need to pass
    a huge mountain to get to the deepest valley. For this kind of problem
    the simulated annealing algorithm is very good. Usualy after a change of
    variables, the combination of parameters with the lowest energy gets picked
    as the state to be in, but if the temperature is high enough sometimes the
    algorithm would also jump to state in which the energy is higher.
    For more details, please look "Simulated annealing" up.


'''

def diag_plot(x_Achse, y_matrix, X_Title, Y_Title, Kurv_Names, PDFTitle, keypos):
    c = graph.graphxy(width=10,height=10,
        x = graph.axis.linear(min = min(x_Achse), max = max(x_Achse),
                          title= X_Title),
        y = graph.axis.linear(min = np.amin(y_matrix),max = np.amax(y_matrix),
                          title= Y_Title),
                      key = graph.key.key(pos=keypos, dist =0.1))
    dd = []
    for i in range(np.shape(y_matrix)[0]):
        y_values = y_matrix[i,:]
        print (len(y_values))
        dd.append(graph.data.values(x = x_Achse, y = y_values ,title = Kurv_Names[i]))

    c.plot(dd,[graph.style.line([color.gradient.Rainbow])])

    c.writePDFfile(PDFTitle)

def diag_plot2(x_Achse, y_Achse, X_Title, Y_Title, Kurv_Name, PDFTitle, keypos):
    c = graph.graphxy(width=10,height=10,
        x = graph.axis.linear(min = min(x_Achse), max = max(x_Achse),
                          title= X_Title),
        y = graph.axis.linear(min = min(y_Achse),max = max(y_Achse),
                          title= Y_Title),
                      key = graph.key.key(pos=keypos, dist =0.1))

    c.plot(graph.data.values(x = x_Achse, y = y_Achse, title  = Kurv_Name),[graph.style.line([color.gradient.Rainbow])])
    c.writePDFfile(PDFTitle)

def diracEigenvalues(n):
    '''
    :param n: Number of current Shell
    '''
    return (2*n + 1)/2


'''
Constraints
To see in TeX-form how the trace-constraint looks like. Look at the
master thesis of nikki kilbertus on page 23 equation 2.55.
The "Constant" in 2.54, and additional algebraic changes after is chosen to be 1.
But can also be varied as one wishes if needed, so it's not fixed to be 1,
but can be changed.

For the integrand i need the different weights rho_1, rho_2 depending on the
shell number N. If N=1, rho_1 is already fixed, because of the constraint.
For N=2 rho_1 or rho_2 can be varied.

Let's say N=3, then i got rho_1 + 3*rho_2 + 6*rho_3 = 1.

So the idea is to just give a random float number to rho_1, rho_2 and rho_3
and then calculate new_constant = 1*5 + 3*7.2 +6*7.5
and then divide everything through new constant.
so we get 1 = (1*5 + 3*7.2 +6*7.5)/new_constant
So rho_1 = 5/new_constant, rho_2 = 7.2/new_constant etc.

For this i at first construct a List with the coefficients of the weights
(Function is called """listmaker"""), which is for N=3, Rho_Koeff_list =
[1,3,6].  The good thing about that is, that the List for bigger N always
contains the smaller one. So N=4, would be Rho_Koeff_list = [1,3,6,10].

In the Variationprocess, for which i will write more in the Variation Function,
i  basically create a List with floats, for example [5, 7.2, 7.5] in the N=3
case (This are the same numbers i used before). This list then is the input for
get_rho_values and in this function, i do exactly what i explained before.  But
be aware that the input list is not completely random and gets changed
according to the variationprocess.
'''

def listmaker(N):
    listim = []
    for n in range(1, N+1):
        listim.append((diracEigenvalues(n)**2 -0.25)/2)

    return listim

def get_rho_values(Factor89, Constant, Rho_Koeffs_List):
    Sum = 0
    for i,rho_koeff in  enumerate(Rho_Koeffs_List):
        Sum += rho_koeff*Factor89[i]
    rho_values = Factor89/Sum
    return rho_values*Constant

def Fitness(Dada, N, x_fitn2):

    print('N', N)
    if Dada ==1:
        return get_Action(t, r, N, Integration_bound, T, *x_fitn2, kappa, False, False, 1)
    if Dada ==2:
        return help_Wirk(N, *x_fitn2)

def which_variant(random_K, random_Rho, random_w):
    '''
    random_K    boolean
    random_Rho  boolean
    random_w    boolean

    It's not always intended to vary all parameters
    at the same time. So this function, will return
    an integeger number, which stands for the variant.
    In x_fitn_func this output is needed to set up
    the initial state.

    This function feeds into
    the initial_state_constructor
    '''
    if random_K == False and random_Rho== False  and random_w== False:
        # This case doesn't make much sense,
        # because of course there is going to be a variation.
        # This case i am using for fixing every paramater initially
        # instead of letting the programm choose a random starting
        # point.
        return 1
    elif random_K == False and random_Rho== False  and random_w==True:
        return 2
    elif random_K == True and random_Rho== False and random_w==False:
        return 3
    elif random_K == False and random_Rho==True  and random_w==False:
        return 4
    elif random_K == False and random_Rho==True  and random_w==True:
        return 5
    elif random_K ==True  and random_Rho== False  and random_w==True:
        return 6
    elif random_K ==True  and random_Rho==True  and random_w==False:
        return 7
    elif random_K ==True  and random_Rho==True  and random_w==True:
        return 8

def Initial_state_constructor(variant, K_List, Rho_List, w_List,
      Rho_Koeffs_List, N):

    '''
    variant
    	integer, input comes from which_variant
    K_List, Rho_List, w_List
    	Lists of length N, with which the initial
        state gets constructed.

    N
    	integer, The Shell number is here needed
        to know which parameter space is being
        searched. Because i could start to probe the
        parameter space from the minimizer, which i
        found on one subspace of the regarded
        parameter space.
        Or just simply i could say, start at this
        point and do not take a random starting point
        (which is default).
        By choosing a starting point and then choosing
        the Shell-Number
        we decide wether we expand into a bigger
        parameterspace or not.
        So if N > len(K_List) then we go into a
        bigger parameter space. If N = len(K_List)
        we just decided to start at some fixed point.
        N < len(K_List) will give an error.


    To set up the initial state, we need to fill up the Lists which are not
    going to get varied.

    This "filling up" can be used to vary with a starting point in mind or
    just keep some half_filled_list fixed.


    '''
    print('erstes N', N)


    Lists_length = len(K_List)
    print(type(Lists_length), Lists_length, type(N))
    if N < Lists_length:
        print("Shell-Number %i < List_length %i. Size-Error. Fix the size of K_List and so on, N can only be equally big or bigger!" %(N, Lists_length) )
        sys.exit()

    half_filled_list = np.zeros((3,N))

    if variant ==1:
        half_filled_list[0, :Lists_length] = K_List
        half_filled_list[1, :Lists_length] = Rho_List
        half_filled_list[2, :Lists_length] = w_List

        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)

    elif variant ==2:
        half_filled_list[0, :Lists_length] = K_List
        half_filled_list[1, :Lists_length] = Rho_List

        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)
    elif variant ==3:
        half_filled_list[1, :Lists_length] = Rho_List
        half_filled_list[2, :Lists_length] = w_List

        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)
    elif variant ==4:
        half_filled_list[0, :Lists_length] = K_List
        half_filled_list[2, :Lists_length] = w_List

        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)
    elif variant ==5:
        half_filled_list[0, :Lists_length] = K_List

        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)
    elif variant ==6:
        half_filled_list[1, :Lists_length] = Rho_List

        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)
    elif variant ==7:
        half_filled_list[2, :Lists_length] = w_List

        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)
    elif variant ==8:
        parameters, ftns_for_x = Gap_Filler(N, Rho_Koeffs_List, half_filled_list)

    return parameters, ftns_for_x



def Gap_Filler(N, Rho_Koeffs_List, half_filled_list):
    if random_K:
        half_filled_list[0]= np.random.random_sample(N)*(K_End - K_Anf)
    if random_Rho:
        rho_randomi = np.random.random_sample(N)
        half_filled_list[1]= get_rho_values(rho_randomi, Constant, Rho_Koeffs_List)
    if random_w:
        half_filled_list[2] = np.random.random_sample(N)*(w_End - w_Anf)
    fitn_wert_x = Fitness(Which_Fitness, N, half_filled_list)

    return half_filled_list, fitn_wert_x

def Variation(N, x_fitn22):
    '''
    N
	integer, Shell-Number
    x_fitn22
	is the array, that contains the parameters of Rho_List,
        pre_w_Listand K_List

    The variation of the state, should not end in a new random state, it
    should be somehow "near" the old state.
    So at the first line of each parameter group, i get a random number
    ranging from -1 to 1 for the weights rho, or -(K_End - K_Anf)/10, + .., .
    The second choice with (K_End - K_Anf)/10 is very arbitrary, but i have
    to start somewhere and it seems to be good.
    The same goes for the frequencies.

    I think everything else is self explanatory.

    The way i programmed the whole programm,
    i only have to change this variation function
    to for exapmple only vary one element of the weights.
    This i do not need right now. But in the future maybe.

    '''

    if var_K:
        randomi2 = (2*np.random.random_sample(N) - 1)*(K_End - K_Anf)/10
        K_randomi5 = np.absolute(x_fitn22[0] + randomi2)
        x_fitn22[0]= list(K_randomi5)

    if var_Rho:
        randomi2 = (2*np.random.random_sample(N) - 1)/10
        rho_randomi6 = np.absolute(x_fitn22[1] + randomi2)/np.array(Rho_Koeffs_List)
        x_fitn22[1] = get_rho_values(rho_randomi6, Constant,
                Rho_Koeffs_List)


    if var_w:
        randomi2 = (2*np.random.random_sample(N) -1)*(w_End - w_Anf)/10
        w_randomi5 = np.absolute(x_fitn22[2] + randomi2)
        x_fitn22[2] = list(w_randomi5)

    return x_fitn22


def Control_Action(N, K_Lte, Rho_Lte, w_Lte):

    '''
    This function serves to test the minimizer.
    As you can see this is a higher dimensional
    parabola. In first order the action is also a
    parabola so if the minimizer won't work here
    it will definitely fail for the action.
    '''
    s = 0
    w_Min = [o for o in range(N)]
    R_Min = [o for o in range(N)]
    K_Min = [o for o in range(N)]

    Min = [w_Min, R_Min, K_Min]
    for ii in range(N):
        if var_K:
            s += (K_Lte[ii] - Min[0][ii])**2
        if var_Rho:
            s += (Rho_Lte[ii] - Min[1][ii])**2
        if var_w:
            s+= (w_Lte[ii] - Min[2][ii])**2
    return s

def K_BoltzmanFinder(Selbst, Rho_Koeffs_List, N):
    '''
    Selbst
        boolean, A value that i choose for the Boltz. const
    Rho_Koeffs_List
        List of floats, I need this for the calculation of
        the Action
    N
        integer, Shell-Number
    var_K
        boolean, Needed to decide whether the Impuls-variables
        should be varied or not
    var_Rho
        boolean, -''-
    var_w
        boolean, -''-
    variant
        integer, it's a number, which is given for every
        kombinatin of var_Rho, var_w, var_w(see
        fucntion which_variant)

    For the simulated annealing process i need a thermal function. This thermal
    function should decrease in time and eventually go to zero.
    For this thermal function i need some sort of exponential decay constant
    (Let's call this the boltzmann constant. You may confuse it with the actual
    boltzmann constant but for now i will call it that),
    and basically that's what i am constructing here.

    There is a problem with this approach. The Boltzmann konstant can be
    calculated to a value which is too high. That is because the process
    just uses somen random number of evaluations of the action and takes
    the average of those to be the boltzmann konstant.

    If the boltzmann konstant is to high, the minimizing algorithm allows
    states to get to higher and higher actions. But the case is, that the
    action in first order is a parabola, one notices that the algorithm may
    just get lost climbing the action mountain. One could just construct a
    better temperatur decay function, but thats all to random. So at some point
    i choose the boltzmann konstant to be 0.001 and sure enough the minimizers
    were more easily found.

    So at some point, i will have to find a better approach for this.
    '''

    K_Boltz = 0
    if Selbst:
        K_Boltz = Selbstval_K_Boltz
    else:
        for _ in range(Mittelgr):
            x_fitn5, fitn_wert_y5 = Initial_state_constructor(variant, K_List,
                    Rho_List, w_List, Rho_Koeffs_List, N)

            K_Boltz += fitn_wert_y5/Mittelgr

    return K_Boltz


def Minimierer(N, first, Rho_Koeffs_List, Candidate_Minimum):

    K_Boltz =K_BoltzmanFinder(True,Rho_Koeffs_List, N)
    kol = open('iterk.txt', 'a')
    #Candidate_Minimum = Initial_State
    print('id(Candidate)',id(Candidate_Minimum[0]))

    fitn_wert_x = Candidate_Minimum[1]#Initial_State[1]
    x_fitn_i = [*Candidate_Minimum[0]]# Initial_State[0]
    iterat = 0
    for m,tt in enumerate(temp):
        for _ in range(4):
            iterat +=1
            new_param_values = Variation (N, x_fitn_i)
            energy_new_param = Fitness(Which_Fitness, N, new_param_values)
            kol.write(str(iterat)+ ' ' + str(new_param_values[0][0]) +' '
                    +str(energy_new_param)+'\n')
            boltzi = boltzmann(fitn_wert_x, energy_new_param, tt, K_Boltz)

            if fitn_wert_x > energy_new_param:
                fitn_wert_x = energy_new_param
                x_fitn_i =  new_param_values
                if Candidate_Minimum[1] > energy_new_param:
                    Candidate_Minimum[0]= new_param_values
                    Candidate_Minimum[1]= energy_new_param
                    print('Cand1', Candidate_Minimum)
            elif (fitn_wert_x < energy_new_param) and (random.random() <= boltzi) :
                fitn_wert_x = energy_new_param
                x_fitn_i =  new_param_values
    kol.close()
    return Candidate_Minimum


if __name__ == "__main__":
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


    #Minimum =[[[0.61402128722400451, 5.4919577589099093], [0.96197069,  0.01267644], [0, 1]], 0.007515003816659801]
    Constant = 1
#    pre2_w_List = eval(pre_w_List)

    pre2_K_List = eval(pre_K_List)

    pre2_Rho_List = eval(pre_Rho_List)

    for SN in range(first,Anzahl_N+1):
        Iter = SN**2                        #Number of temperatur iterations
        hilfsarray_fuer_temp = np.linspace(0.01,5,Iter)
        Amplitude = 0.1                     #Amplitude of tempearatur oszillation
                                            #on the exponentially decreasing
                                            #temperatur
        freq = np.pi                        #Frequenz for oscillation
        halb = 0.001                        #exponential decay constant
        temp = temperatur(Iter, hilfsarray_fuer_temp, halb, freq, Amplitude)
        Factor = np.random.random_sample(SN)
        Rho_Koeffs_List = listmaker(SN)

        pre2_w_List = eval(pre_w_List) # I put this list assignment here,
                                    #    should be right. Not entirely sure,
                                    #    but i want to let the code run for the
                                    #    next hours and i am in a rush.

        x_fitn11, fitn_wert_y11 = Initial_state_constructor(variant,
           pre2_K_List, pre2_Rho_List, pre2_w_List, Rho_Koeffs_List, SN)


        print('Heyhooo', x_fitn11)
        Initial_State = [x_fitn11, fitn_wert_y11]
        print('Initial_State', Initial_State)
        Minimum = Minimierer(SN, first, Rho_Koeffs_List, Initial_State)

        #pre2_w_List = [*Minimum[0][2],0]
        print('pre2_w_list', pre2_w_List)
        pre2_K_List = [*Minimum[0][0],0]
        pre2_Rho_List = [*Minimum[0][1],0]

        gg = open('Minimum7.txt', 'a')
        gg.write('Minimum fuer N = %d'%(SN) + str(Minimum)+'\n')
        gg.close()
