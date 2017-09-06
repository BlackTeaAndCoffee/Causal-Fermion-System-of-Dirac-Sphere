from sympy.printing import ccode
from scipy import integrate
from symengine import I
from sympy import *
from configfunktion import configfunktion
import configparser
import sympy as sy
import symengine as si
import numpy as np
import random
import time
import os
import sys
import Rho_data
import scipy.misc
import ctypes
from SymEngineFast import *
from SimTest import *

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
     return (2*n + 1)/2


'''
Constraints
'''
def traceConstraint(N, Constant, rho):
    #dann benutze rho[0] usw.
    addd = Constant
    for n in range(1, N):
        addd -= rho[n-1]*(diracEigenvalues(n)**2 - 0.25)

    return addd


'''
Constraints
'''

def subs_coeffs(N, Constant, rho):
    subs  = traceConstraint(NN, Constant, rho)
    liste = []
    if N == 1:
        return [1]
    else:
        for ii in range(N - 1):
            liste.append(Poly(subs).coeff_monomial(rho[ii]))
    return liste

def listmaker(N, Constant):
    listim = []
    for ni in range(1, N+1):
        listim.append((diracEigenvalues(ni)**2 -0.25)/2)

    return listim

def get_rho_values2(N, Factor89,for_rho, Const, Rho_Koeffs_List, SameValues=True):
    Sum = 0
    for i in range(N):
        Sum += Rho_Koeffs_List[i]*Factor89[i]
    rho_values = Factor89/Sum
    return rho_values

def Fitness(Dada, N, x_fitn2):
    if Dada ==1:
        return get_Action(t, r, N, Intgrenze, T, *x_fitn2, kappa, False, False, 1)
    if Dada ==2:
        return help_Wirk(N, *x_fitn2)
def Zeilensparer(Rho_Koeffs_List, Minimum2, N,first, x_fitn2, var_K, var_Rho,var_w,
        variant, SartWithGivenMinima):

    if SartWithGivenMinima and N != 1 :  #Diese SartWithGivenMinima ist dafuer da, dass ich beim testen mit der
                  #HilfsAction ohne auf die Falle weiter unten einzugehen
                  #fortfahren kann.
        #print(x_fitn2)
        if var_K:
            x_fitn2[0]= [*Minimum2[0][0],0]
        if var_Rho:
            x_fitn2[1]= [*Minimum2[0][1],0]
        if var_w:
            x_fitn2[2]= [*Minimum2[0][2],0]
        #print('lala', x_fitn2)

        fitn_wert_y2 = Fitness(Fitness_Nr, N, x_fitn2)#help_Wirk(N, *x_fitn2)#get_Action(t, r, N, Intgrenze, T, *x_fitn2, kappa, False, False, 1)
    else:
        if var_K:
            x_fitn2[0]= np.random.random_sample(N)*(K_End - K_Anf)
        if var_Rho:
            rho_randomi = np.random.random_sample(N)
            x_fitn2[1]= get_rho_values2(N,rho_randomi,0, Constant, Rho_Koeffs_List, SameValues  =False)
        if var_w:
            x_fitn2[2] = np.random.random_sample(N)*(w_End - w_Anf)
        fitn_wert_y2 = Fitness(Fitness_Nr, N, x_fitn2)#help_Wirk(N, *x_fitn2)#get_Action(t, r, N, Intgrenze, T, *x_fitn2, kappa, False, False, 1)
    return x_fitn2, fitn_wert_y2

def x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste):
    if variant ==1:
        x_fitn5 = [K_Liste, Rho_Liste, w_Liste]
    elif variant ==2:
        x_fitn5 = [K_Liste, Rho_Liste, []]
    elif variant ==3:
        x_fitn5 = [[], Rho_Liste, w_Liste]
    elif variant ==4:
        x_fitn5 = [K_Liste, [], w_Liste]
    elif variant ==5:
        x_fitn5 = [K_Liste, [], []]
    elif variant ==6:
        x_fitn5 = [[], Rho_Liste, []]
    elif variant ==7:
        x_fitn5 = [[],[],w_Liste]
    elif variant ==8:
        x_fitn5 = [[],[],[]]

    return x_fitn5

def Variierer(N, x_fitn22, var_rho = True, var_K = True, var_w = True):

#       K_randomi = np.random.random_sample(N)*(K_End - K_Anf)
#       w_randomi = np.random.random_sample(N)*(w_End - w_Anf)

#       K_randomi = np.absolute(x_y_min[0][0] + (2*np.random.random_sample(N)
#           - 1)*(K_End - K_Anf))
#       print('K_randomi = ', K_randomi)
#       w_randomi = np.absolute(x_y_min[0][2] + (2*np.random.random_sample(N)
#           - 1)*(w_End - w_Anf))
#       rho_randomi = np.random.random_sample(N)
#       for_rho = np.random.randint(0,N)

    if var_K:
        randomi2 = (2*np.random.random_sample(N) - 1)*(K_End - K_Anf)/10
        K_randomi5 = np.absolute(x_fitn22[0] + randomi2)
        x_fitn22[0]= list(K_randomi5)
#   if var_rho:
#       randomi2 = (2*np.random.random() - 1)/10
#       rho_randomi6 = np.absolute(rho_randomi5 + randomi2)
#       while rho_randomi6 > 1:
#           randomi2 = (2*np.random.random() - 1)/10
#           rho_randomi6 = np.absolute(rho_randomi5 + randomi2)
#       x_fitn[1] = get_rho_values2(N,rho_randomi6,for_rho,Constant,Rho_Koeffs_List, SameValues  =False)
#       print('rho_randomi5 =', rho_randomi5)
    if var_rho:
        randomi2 = (2*np.random.random_sample(N) - 1)/10
        rho_randomi6 = np.absolute(x_fitn22[1] + randomi2)/np.array(Rho_Koeffs_List)
        x_fitn22[1] = get_rho_values2(N,rho_randomi6,for_rho,Constant,Rho_Koeffs_List, SameValues  =False)


    if var_w:
        randomi2 = (2*np.random.random_sample(N) -1)*(w_End - w_Anf)/10
        w_randomi5 = np.absolute(x_fitn22[2] + randomi2)
        x_fitn22[2] = list(w_randomi5)

    return x_fitn22

def which_variant(var_K, var_Rho, var_w):
    if var_K == False and var_Rho== False  and var_w== False:
        return 1
    elif var_K == False and var_Rho== False  and var_w==True:
        return 2
    elif var_K == True and var_Rho== False and var_w==False:
        return 3
    elif var_K == False and var_Rho==True  and var_w==False:
        return 4
    elif var_K == False and var_Rho==True  and var_w==True:
        return 5
    elif var_K ==True  and var_Rho== False  and var_w==True:
        return 6
    elif var_K ==True  and var_Rho==True  and var_w==False:
        return 7
    elif var_K ==True  and var_Rho==True  and var_w==True:
        return 8

def help_Wirk(N, K_Lte, Rho_Lte, w_Lte):
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

def K_BoltzmanFinder(Selbst, Rho_Koeffs_List, N, first, var_K, var_Rho,
        var_w, variant, StartWithGivenMinima):
    K_Boltz = 0
    if Selbst:
        K_Boltz = Selbstval_K_Boltz
    else:
        for _ in range(Mittelgr):
            x_fitn5 = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
            x_fitn5, fitn_wert_y5 = Zeilensparer(Rho_Koeffs_List, Minimum, N,
                    first, x_fitn5, var_K, var_Rho, var_w, variant, False)

            K_Boltz += fitn_wert_y5/Mittelgr

    return K_Boltz


def Minimierer(N, first, Rho_Koeffs_List, var_K, var_Rho, var_w,
         variant, StartWithGivenMinima, x_y_min):

    K_Boltz =K_BoltzmanFinder(True,Rho_Koeffs_List, N, first, var_K, var_Rho,
            var_w, variant, StartWithGivenMinima)
    kol = open('iterk.txt', 'a')
    x_y_min = Initial_State
    fitn_wert_x = Initial_State[1]
    x_fitn_i = Initial_State[0]
    K_Iter_List = [x_fitn_i[0][0]]
    m_list = [0,1,2,3]

    iterat = 0
    for m,tt in enumerate(temp):
#       if m == m_list[m]: #Diese If-Abfrage existiert weil, ich zuerst nur in
#                          #der Naehe des Minums variieren moechte.
#           K_randomi = [*Minimum[0],0]
#           w_randomi = [*Minimum[1],0]
#           rho_randomi = [*Minimum[2],0]
#       else:
#           K_randomi = np.random.random_sample(N)*(K_End - K_Anf)
#           w_randomi = np.random.random_sample(N)*(w_End - w_Anf)

#       K_randomi = np.absolute(x_y_min[0][0] + (2*np.random.random_sample(N)
#           - 1)*(K_End - K_Anf))
#       print('K_randomi = ', K_randomi)
#       w_randomi = np.absolute(x_y_min[0][2] + (2*np.random.random_sample(N)
#           - 1)*(w_End - w_Anf))
#           rho_randomi = np.random.random_sample(N)
#       for_rho = np.random.randint(0,N)
        for _ in range(4):
            iterat +=1
            x_fitn33 = Variierer (N, x_fitn_i, var_Rho, var_K, var_w)
            fitn_wert_y = Fitness(Fitness_Nr, N, x_fitn33)#help_Wirk(N, *x_fitn) #get_Action(t, r, N, Intgrenze, T,
                    #*x_fitn, kappa, False, False, 1)
#           K_randomi = x_fitn[0]
#           w_randomi  = x_fitn[2]
            kol.write(str(iterat)+ ' ' + str(x_fitn33[0][0]) +' '
                    +str(fitn_wert_y)+'\n')
            print('VArs and y_val', x_fitn33, fitn_wert_y)
            boltzi = boltzmann(fitn_wert_x, fitn_wert_y, tt, K_Boltz)
            if fitn_wert_x > fitn_wert_y:
                fitn_wert_x=fitn_wert_y
                if x_y_min[1] > fitn_wert_y:
                    x_y_min[0]= x_fitn33
                    x_y_min[1]= fitn_wert_y

            elif (fitn_wert_x < fitn_wert_y) and (random.random() <= boltzi) :
                print(fitn_wert_x, fitn_wert_y, 'jump')
                fitn_wert_x = fitn_wert_y
    kol.close()
    return x_y_min

if __name__ == "__main__":
    '''
    Variablen für das kausale System werden nachfolgend deklariert
    '''

    var_K, var_Rho, var_w = configfunktion('Vary_Parameters_bool')
    K_Anf, K_End, K_List = configfunktion('Impuls')
    w_Anf, w_End, w_List = configfunktion('Frequenz')
    Constant, kappa, Rho_List = configfunktion('Constraints')
    Anzahl_N, first = configfunktion('System_sizes')
    StartWithGivenMinima, Fitness_Nr = configfunktion('Test')
    variant = which_variant(var_K, var_Rho, var_w)

    print(StartWithGivenMinima)

    T = 2*np.pi
    r = si.symarray('r', 1)
    t = si.symarray('t', 1)

#   K_Anf = 0
#   K_End = 25
#   K_Liste = list(np.linspace(K_Anf, K_End, 10))


#   kappa = 0.0001
#   kappa_Anzahl = 1


#   #Rho_Liste = [0.8,0.07]#Rho_data.get_rho_values2(N,Factor, Constant, Rho_Koeffs_List, SameValues  = True)

    x_Anf = 0
    x_End = np.pi

    Intgrenze = [x_Anf, x_End]
    Wirk = []
    Selbstval_K_Boltz =0.001
    Mittelgr = 4
#   w_Anf = 0
#   w_End = 10
    for_rho = 1

    '''
    Ab hier werden die fuer die Minimierung benötigten Variablen definiert
    '''
    '''
    Ich passe das Minimieren zuerst fuer den eindim. Fall an, danach
    passe ich es für mehrere dim. an.
    '''
#   Anzahl_N = 4
#   Minimum  = [[9.7903003749135973, 1.0428185324876262], np.array([ 0.24596702,
#       0.25134433]), [0, 1]] #kappa = 10-4
    Minimum =[[[0.61402128722400451, 5.4919577589099093], [0.96197069,  0.01267644], [0, 1]], 0.007515003816659801] # kappa = 0.5

    #   first = 3
    for SN in range(first,Anzahl_N+1):
        Iter = 20*SN                                 #Die Anzahl der temperaturiterationen
        hilfsarray_fuer_temp = np.linspace(0.01,5,Iter)
        Amplitude = 0.1                                #Amplitude der Fluktuation
                                                   #der Temperatur
        freq = np.pi                                    #Die Frequenz mit der die temp variieren soll
        halb = 0.001                                      #Der Halbwertswert fuer die
        temp = temperatur(Iter, hilfsarray_fuer_temp, halb, freq, Amplitude)


        w_Liste = eval(w_List)
        print('w_Liste', type(w_Liste))
        Rho_Koeffs_List = listmaker(SN, Constant)

        #rho_randomi = j*0.1
        #Rho_Liste =get_rho_values2(j,rho_randomi, Constant, liste, SameValues  =False)
        #print('Rho_Liste', Rho_Liste)
        K_Liste = eval(K_List)
        Constant = 1
        Factor = np.random.random_sample(SN)
        print('Rho_Koeffs_List', Rho_Koeffs_List)
        if var_Rho == True:
            Rho_Liste = get_rho_values2(SN, Factor, 1, Constant, Rho_Koeffs_List, SameValues  = False)
            print('Rho', Rho_Liste)

        else:
            Rho_Liste = eval(Rho_List)

        print('Rho_Liste = ', Rho_Liste)

        x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
        print('Heyhooo', x_fitn)
        x_fitn11, fitn_wert_y11 = Zeilensparer(Rho_Koeffs_List, Minimum, SN, first, x_fitn,
            var_K, var_Rho,var_w, variant,StartWithGivenMinima)

        print('Heyhooo', x_fitn11)
        Initial_State = [x_fitn11, fitn_wert_y11]
        print('Initial_State', Initial_State)
        Minimum = Minimierer(SN, first, Rho_Koeffs_List,var_K,
                var_Rho, var_w, variant, StartWithGivenMinima,
                Initial_State)

        gg = open('Minimum7.txt', 'a')
        gg.write('Minimum fuer N = %d'%(SN) + str(Minimum)+'\n')
        gg.close()
