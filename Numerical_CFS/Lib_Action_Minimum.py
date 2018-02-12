from sympy.printing import ccode
from scipy import integrate
from symengine import I
from sympy import *
from Numerical_CFS.configfunktion import configfunktion
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
#from .LibForSimulAnnealing import *

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

    * There is another method called bayesian optimization. I haven't tried
      this one out yet. But i will definitely do so soon.

In the simulated annealing algorithm we basically begin with some initial
state. A state is just some combinations of the various paramters in
K_List, Rho_List and w_List.

    * The initial state should either be given by the user or the programm should
      pick some state.

    * The user would give a state, if he or she want's to give the algorithm
      a starting point. This we will do by each increase of the Shell-Number
      N.

    * A state always gets an Energyvalue (float number) and the state with the
      lowest Energy is the desired one. Now we can imagine, that the energy
      landscape has a lot of mountains and valleys and sometime you need to pass
      a huge mountain to get to the deepest valley. For this kind of problem
      the simulated annealing algorithm is very good. Usualy after a change of
      variables, the combination of parameters with the lowest energy gets picked
      as the state to be in, but if the temperature is high enough sometimes the


'''

def diag_plot(x_Axis, y_matrix, X_Title, Y_Title, Curve_Names, PDF_Name, keypos):
    '''
    I use `PyX <http://pyx.sourceforge.net>`_ for plotting.

    :param x_Axis: A List with the values of the x-Axis
    :param y_Matrix: A (Number of Curves)x(Number of elements in x_Axis) Matrix.
    :param X_Title: Title of x-Axis.
    :param Y_Title: Title of y-Axis.
    :param Curve_Names: You can name each and every curve.
    :param PDF_Name: The name of the PDF.
    :param keypos: Position of the list of Curvenames in the PDF. You can choose\
    "t" for top, "l" for left and "r" for right, "b" for bottom, "c" for center \
    and "m" for middle. You can also use a combination. For example "tr".
    :type x_Axis: A numpy.array, the length of is up to you.
    :type y_Matrix: A numpy.array.
    :type X_Title: A String.
    :type Y_Title: A String.
    :type Curve_Names: List of Strings.
    :type PDF_Name: A String.
    :type keypos: A String.
    :return: The programm creates a file with "PDF-Title" as it's name in the same\
    directory you have run this programm in.
    '''


    c = graph.graphxy(width=10,height=10,
        x = graph.axis.linear(min = min(x_Axis), max = max(x_Axis),
                          title= X_Title),
        y = graph.axis.linear(min = np.amin(y_matrix),max = np.amax(y_matrix),
                          title= Y_Title),
                      key = graph.key.key(pos=keypos, dist =0.1))
    dd = []
    for i in range(np.shape(y_matrix)[0]):
        y_values = y_matrix[i,:]
        print (len(y_values))
        dd.append(graph.data.values(x = x_Axis, y = y_values ,title = Curve_Names[i]))

    c.plot(dd,[graph.style.line([color.gradient.Rainbow])])

    c.writePDFfile(PDF_Name)



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

In the type_1_variationprocess, for which i will write more in the type_1_variation Function,
i  basically create a List with floats, for example [5, 7.2, 7.5] in the N=3
case (This are the same numbers i used before). This list then is the input for
get_rho_values and in this function, i do exactly what i explained before.  But
be aware that the input list is not completely random and gets changed
according to the variationprocess.
'''
class Initial_System_Params():
    def __init__(self, random_K, random_Rho, random_w, K_List, Rho_List, w_List, kappa):
        self.random_K = random_K 
        self.random_Rho = random_Rho
        self.random_w = random_w
        self.K_List = K_List
        self.w_List = w_List
        self.Rho_List = Rho_List
        print(w_List)
        self.N = len(w_List) #w_List, K_List and Rho_List have all the same length.
        print('N =',self.N)
        self.variant = self.which_variant() 
        self.kappa = kappa
    def which_variant(self):
        '''
        :param random_K: This paramter gets its default value from settings.cfg.
        :param random_Rho: This paramter gets its default value from settings.cfg.
        :param random_w: This paramter gets its default value from settings.cfg.
        :type random_K:  boolean
        :type random_Rho: boolean
        :type random_K: boolean

        It's not always intended to vary all parameters.\
        Some just fixed values.\
        So this function, will return\
        an integer number, which stands for the variant.\
        In x_fitn_func this output is needed to set up\
        the initial state.\

        This function feeds into\
        the initial_state_constructor.
        '''
        if self.random_K == False and self.random_Rho== False  and self.random_w== False:
            # This case doesn't make much sense,
            # because of course there is going to be a variation.
            # This case i am using for fixing every paramater initially
            # instead of letting the programm choose a random starting
            # point.
            return 1
        elif self.random_K == False and self.random_Rho== False  and self.random_w==True:
            return 2
        elif self.random_K == True and self.random_Rho== False and self.random_w==False:
            return 3
        elif self.random_K == False and self.random_Rho==True  and self.random_w==False:
            return 4
        elif self.random_K == False and self.random_Rho==True  and self.random_w==True:
            return 5
        elif self.random_K == True and self.random_Rho== False  and self.random_w==True:
            return 6
        elif self.random_K == True and self.random_Rho==True  and self.random_w==False:
            return 7
        elif self.random_K == True and self.random_Rho==True  and self.random_w==True:
            return 8

    def Initial_Params_Constructor(self):

        '''
        :param variant: input comes from which_variant
        :type variant: Integer
        :param K_List: These impulse-variables are one part of the set of variables,\
        which later are needed to be choosen new every time before calculating the action.
        :type K_List: List of length N.
        :param w_List: These frequence-variables are one part of the set of variables,\
        which later are needed to be choosen new every time before calculating the action.
        :type w_List: List of length N.
        :param Rho_List: These Weight-variables are one part of the set of variables,\
        which later are needed to be choosen new every time before calculating the action.
        :type Rho_List: List of length N.
        :param Rho_Koeffs_List(SN): This List is needed for the calculation of Rho_List.\
        The weights are a special case, because they need to obey certain Constraints.
        :type Rho_Koeffs_List(SN): List of length N.
        :param N: The Shell number is here needed\
            to know which parameter space is being\
            searched. Because i could start to probe the\
            parameter space from the minimizer, which i\
            found on one subspace of the regarded\
            parameter space.\
            Or just simply i could say, start at this\
            point and do not take a random starting point\
            (which is default).\
            By choosing a starting point and then choosing\
            the Shell-Number\
            we decide wether we expand into a bigger\
            parameterspace or not.\
            So if N > len(K_List) then we go into a\
            bigger parameter space. If N = len(K_List)\
            we just decided to start at some fixed point.\
            N < len(K_List) will give an error.\

        :type N: Integer.

        To set up the initial state, we need to fill up the Lists which are not\
        going to get varied.\

        This "filling up" can be used to set the starting point or\
        just keep some half_filled_list as it is. Later on the code \
        will fill the gaps.


        '''
        print('erstes N', self.N)


        Lists_length = self.N
        print(type(Lists_length), Lists_length, type(self.N))
        if self.N < Lists_length:
            raise ValueError("Shell-Number %i < List_length %i. Size-Error. Fix the size of K_List and so on, N can only be equally big or bigger!" %(N, Lists_length))
            sys.exit()

        half_filled_list = np.zeros((3,self.N))

        if self.variant ==1:
            half_filled_list[0, :Lists_length] = self.K_List
            half_filled_list[1, :Lists_length] = self.Rho_List
            half_filled_list[2, :Lists_length] = self.w_List

            parameters = self.Gap_Filler(half_filled_list)

        elif self.variant ==2:
            half_filled_list[0, :Lists_length] = self.K_List
            half_filled_list[1, :Lists_length] = self.Rho_List

            parameters = self.Gap_Filler(half_filled_list)
        elif self.variant ==3:
            half_filled_list[1, :Lists_length] = self.Rho_List
            half_filled_list[2, :Lists_length] = self.w_List

            parameters = self.Gap_Filler(half_filled_list)
        elif self.variant ==4:
            half_filled_list[0, :Lists_length] = self.K_List
            half_filled_list[2, :Lists_length] = self.w_List

            parameters = self.Gap_Filler( half_filled_list)
        elif self.variant ==5:
            half_filled_list[0, :Lists_length] = self.K_List

            parameters = self.Gap_Filler(half_filled_list)
        elif self.variant ==6:
            half_filled_list[1, :Lists_length] = self.Rho_List

            parameters = self.Gap_Filler(half_filled_list)
        elif self.variant ==7:
            half_filled_list[2, :Lists_length] = self.w_List

            parameters = self.Gap_Filler(half_filled_list)
        elif self.variant ==8:
            parameters = self.Gap_Filler(half_filled_list)

        return [*parameters , self.kappa]

    def Gap_Filler(self, half_filled_list):
        if self.random_K:
            half_filled_list[0]= np.random.random_sample(self.N)*(K_End - K_Anf)
        if self.random_Rho:
            rho_randomi = np.random.random_sample(self.N)
            half_filled_list[1]= Rho_values(rho_randomi)
        if self.random_w:
            half_filled_list[2] = np.random.random_sample(self.N)*(w_End - w_Anf)
        return half_filled_list

class Variation_of_Parameters():
    def __init__(self, var_K, var_Rho, var_w,  delta_K, delta_Rho, delta_w, Rho_Values):
        self.var_K = var_K 
        self.var_Rho = var_Rho
        self.var_w = var_w
        self.delta_K = delta_K
        self.delta_w = delta_w
        self.delta_Rho = delta_Rho
        self.Rho_Values = Rho_Values 
    def __call__(self, Input_List):
        '''
        :param N: Integer.
        :type N: Shell-Number.
        :param self.Input_List: Is the array, that contains the parameters of Rho_List,\ 
                         pre_w_Listand K_List.
        :type  self.Input_List: List of 3 Lists.
        :param var_K: This is a global variable. We need them in this\
        function. You need to decide whether the Impuls-variables should be varied or not.
        :type var_K: boolean.
        :param var_Rho: This is a global variable. We need them in this\
        function. You need to decide whether the Weight-variables should be varied or not.
        :type var_Rho: boolean.
        :param var_w: This is a global variable. We need them in this\
        function. You need to decide whether the frequence-variables should be varied or not.
        :type var_w: boolean.
        :param variant: It's a number, which is given for every kombinatin of var_Rho, 
                        var_w, var_w(see function which_variant).\

        The variation of the state, should not end in a new random state, it\
        should be somehow "near" the old state.\
        So at the first line of each parameter group, i get a random number\
        ranging from -1 to 1 for the weights rho, or -(K_End - K_Anf)/10, + .., .\
        The second choice with (K_End - K_Anf)/10 is very arbitrary, but i have\
        to start somewhere and it seems to be good.\
        The same goes for the frequencies. \
        The way i programmed the whole programm,\
        i only have to change this variation function\
        to for exapmple only vary one element of the weights.\
        This i do not need right now. But in the future maybe.\
        '''
        N = len(Input_List[0])
        if self.var_K:
            randomi2 = (2*np.random.random_sample(N) - 1)*self.delta_K
            K_randomi5 = np.absolute(Input_List[0] + randomi2)
            Input_List[0]= list(K_randomi5)
            print('InpoutLiust', Input_List)
        if self.var_Rho:
            randomi2 = (2*np.random.random_sample(N) - 1)*self.delta_Rho
            rho_randomi6 = np.absolute(Input_List[1] + randomi2)
            Input_List[1] = self.Rho_Values(rho_randomi6)

        if self.var_w:
            randomi2 = (2*np.random.random_sample(N) -1)*self.delta_w
            w_randomi5 = np.absolute(Input_List[2] + randomi2)
            Input_List[2] = list(w_randomi5)
 
        return Input_List

class Rho_Class:
    def __init__(self, N, Constant):
        self.N = N
        self.Constant = Constant
        self.Rho_List = self.Rho_Koeffs_List()
    def __call__(self, Incoming_List):
        Sum = 0
        for i,rho_koeff in  enumerate(self.Rho_List):
            Sum += rho_koeff*Incoming_List[i]
        rho_values = Incoming_List/Sum
        return rho_values*self.Constant
    
    def Rho_Koeffs_List(self):
        listim = []
        for n in range(1, self.N + 1):
            listim.append((self.diracEigenvalues(n)**2 -0.25)/2)

        return listim

    def SameRhoValues(self):
        IdenticalRhos = np.ones(self.N)*self.Constant
        summe = 0
        for i in range(1, self.N + 1):
            summe += ((self.diracEigenvalues(i))**2 - 0.25)/2
        
        return list(IdenticalRhos*(1/summe))

    def diracEigenvalues(self, n):
        '''
        :param n: Number of current Shell
        '''
        return (2*n + 1)/2



class Simulated_Annealing():
    def __init__(self, BaseArrayForTemp, Boltzmann_Constant, 
                     decay_constant, freq, Amplitude, vary, Fitness):
        self.Boltzmann_Constant = Boltzmann_Constant
        self.BaseArrayForTemp = BaseArrayForTemp
        self.decay_constant = decay_constant
        self.freq = freq
        self.Amplitude = Amplitude
        self.vary = vary
        self.Fitness = Fitness
    def boltzmann(self, f_x, f_y, temp):
        '''
        This function will return for an Energydifferenz = f_x - f_y, the boltzmann\
        probability.
        '''
        
        print('f_x,f_y', f_x, f_y)
        diff = abs(f_y - f_x)
        print('diff', diff)
        return np.exp(- diff/(temp*self.Boltzmann_Constant))

    def temperatur_Function(self, temp_iter):
        '''
        This function i need for nonlinear Temperatur-curves.
        :return: Temperatur that is nonlinear.. . this feeds into temperatur, which\ 
        produces the list of temperatures.
        '''
        return np.exp(- self.decay_constant*temp_iter**2)*(self.Amplitude*np.cos(self.freq*temp_iter*0.5) + 1)

    def temperatur(self):
        '''
        The Temperatur which shall be put in the boltzmann-function
        gets put together here.
        '''
        temperatur_list =[]
        for heat_value in self.BaseArrayForTemp:
            temperatur_list.append(self.temperatur_Function(heat_value))
        return temperatur_list

    def Minimierer(self, Candidate_Minimum):
        N = len(Candidate_Minimum[0][1])
        kol = open('iterk.txt', 'a')
        #Candidate_Minimum = Initial_State
        print('id(Candidate)',id(Candidate_Minimum[0]))

        fitn_wert_x = Candidate_Minimum[1]#Initial_State[1]
        x_fitn_i = [*Candidate_Minimum[0]]# Initial_State[0]
        iterat = 0
        temp = self.temperatur()
        for m,tt in enumerate(temp):
            for _ in range(4):
                iterat +=1
                new_param_values = self.vary(x_fitn_i)
                
                self.Fitness.K_Liste = new_param_values[0] 
                self.Fitness.Rho_Liste = new_param_values[1]
                self.Fitness.w_Liste = new_param_values[2]
                self.Fitness.kappa = new_param_values[3]
                        
                energy_new_param = self.Fitness.get_Action()
                kol.write(str(iterat)+ ' ' + str(new_param_values[0][0]) +' '
                        +str(energy_new_param)+'\n')
                boltzi = self.boltzmann(fitn_wert_x, energy_new_param, tt)

                if fitn_wert_x > energy_new_param:
                    fitn_wert_x = energy_new_param
                    x_fitn_i =  new_param_values
                    if Candidate_Minimum[1] > energy_new_param:
                        Candidate_Minimum[0]= [*new_param_values]
                        Candidate_Minimum[1]= energy_new_param
                        print('Cand1', Candidate_Minimum)
                elif (fitn_wert_x < energy_new_param) and (random.random() <= boltzi) :
                    fitn_wert_x = energy_new_param
                    x_fitn_i =  new_param_values
        kol.close()
        print('Candidate_Minimum_adsfadfa', Candidate_Minimum)
        return Candidate_Minimum

   
def MainProg():
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


    #Minimum =[[[0.61402128722400451, 5.4919577589099093], [0.96197069,  0.01267644], [0, 1]], 0.007515003816659801]
    Constant = 1
#    pre2_w_List = eval(pre_w_List)

    pre2_K_List = eval(pre_K_List)

    pre2_Rho_List = eval(pre_Rho_List)

    
    for SN in range(first, Anzahl_N+1):
        Iter = 10 + SN**2                        #Number of temperatur iterations
        BaseArrayForTemp = np.linspace(0.01,5,Iter)
        Amplitude = 0.1                     #Amplitude of tempearatur oszillation
                                            #on the exponentially decreasing
                                            #temperatur
        freq = np.pi                        #Frequenz for oscillation
        decay_constant = 0.001                        #exponential decay constant
 
        Rho_Values = Rho_Class(SN, Constant)
        vary = Variation_of_Parameters(var_K, var_Rho, var_w, delta_K, delta_Rho, delta_w, Rho_Values)
        pre2_w_List = eval(pre_w_List) # I put this list assignment here,
                                       #I set the list in settings.cfs, and it's like 
                                       #[i for i in range(SN)]. So i need SN.
        print('w_List', pre2_w_List)
        Sys_Params= Initial_System_Params(random_K, random_Rho, random_w, 
            pre2_K_List, pre2_Rho_List, pre2_w_List, kappa)
    
        variant = Sys_Params.variant

        print(StartWithGivenMinima, "variant = ", variant)

 

        System_Parameters= Sys_Params.Initial_Params_Constructor()
        
        print('System_Parameters =', System_Parameters)
        CFS_Action = C_F_S(SN, Integration_bound, T, System_Parameters, Schwartzfunktion = True, 
        Comp_String = False, Integration_Type = 1, Test_Action = False)
        Minimum_Finder = Simulated_Annealing(BaseArrayForTemp, Boltzmann_Constant, 
                            decay_constant, freq, Amplitude, vary, CFS_Action)


        fitn_wert_y11 = CFS_Action.get_Action()

    
        Initial_State = [System_Parameters, fitn_wert_y11]
        print('Initial_State', Initial_State)
  
        Minimum = Minimum_Finder.Minimierer(Initial_State)

        #pre2_w_List = [*Minimum[0][2],0]
        print('pre2_w_list', pre2_w_List)
        pre2_K_List = [*Minimum[0][0],0]
        pre2_Rho_List = [*Minimum[0][1],0]

        gg = open('Minimum7.txt', 'a')
        gg.write('Minimum fuer N = %d'%(SN) + str(Minimum)+'\n')
        gg.close()


if __name__ == "__main__":
    MainProg()
