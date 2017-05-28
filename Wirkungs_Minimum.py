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
import Rho_data
import scipy.misc
import ctypes
from SymEngineFast import *
from SimTest import *

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
    add = Constant
    for n in range(1, N):
        add -= rho[n-1]*(diracEigenvalues(n)**2 - 0.25)

    return add


'''
Constraints
'''

def subs_coeffs(NN, Constant, rho):
    subs  = traceConstraint(NN, Constant, rho)
    liste = []
    if NN == 1:
        return [1]
    else:
        for i in range(NN - 1):
            liste.append(Poly(subs).coeff_monomial(rho[i]))
    return liste

def get_rho_values(NN, Factor, for_rho, Constant, liste_in, SameValues = True):
    print('llist', liste_in)
    mix_list = np.arange(NN-1)
    print('asdfasdf',mix_list)
    np.random.shuffle(mix_list)
    print('asfasdf',mix_list)

    if SameValues:
        s = 0
        for i in range(1, NN+1):
            s += (diracEigenvalues(i)**2 - 0.25)
        rho_values = [Constant/s for i in range(NN)]
    else:
        wert = Constant
        rho_values = np.zeros(NN)
        if NN ==1:
            return [2]
        elif NN ==2:
            subb = Factor*wert
            rho_values[0] = subb/liste_in[0]
            print('rho_val', rho_values[0], liste_in)
            wert = wert - subb
            if wert == 0:
                return rho_values

            rho_values[NN-1] = 2*wert/(diracEigenvalues(NN)**2-0.25)
        else:
            print('wert', wert)
            subb2 = Factor*wert
            rho_values[for_rho] = Factor*wert/liste_in[for_rho]
            print('rho_vals', rho_values)
            wert = wert - rho_values[for_rho]*liste_in[for_rho]


            for jl,listval  in enumerate(mix_list):
                print('wert', wert)
                if wert ==0:
                    return rho_values
                elif listval == for_rho:
                    continue
                else:
                    a = random.random()
                    subb3 = a*wert
                    rho_values[listval] = a*wert/liste_in[listval]    # hier steht ein Minus,weil die
                                                            # Koeffizienten negativ sind,
                                                           # finden aber als positive
                                                           # Zahlen Verwendung.
                    print('zahl', wert)
                    wert = wert - rho_values[listval]*liste_in[listval]
                    print('zahl', wert)
                    if wert == 0:
                        return rho_values

            rho_values[NN-1] = 2*wert/(diracEigenvalues(NN)**2 - 0.25)
    print('rho_values, liste_in :', rho_values, liste_in)
    s = 0

    for ll in range(len(rho_values)-1):
        s+= liste_in[ll]*rho_values[ll]
        print(s)
    s+= rho_values[-1]*(diracEigenvalues(NN)**2 - 0.25)/2
    print('Summe',s)
    return rho_values

def Zeilensparer(liste2, Minimum, N,first, x_fitn, var_K, var_Rho,var_w, variant):
    print('x_fitn', x_fitn, liste)

    if first ==N:
        if var_K:
            x_fitn[0]= np.random.random_sample(N)*(K_End - K_Anf)
        if var_Rho:
            rho_randomi = np.random.random()
            x_fitn[1]= get_rho_values(N,rho_randomi,0, Constant, liste2, SameValues  =False)
        if var_w:
            x_fitn[2] = np.random.random_sample(N)*(w_End - w_Anf)
        fitn_wert_y =  get_Wirkung(t, r, N, Intgrenze, T, *x_fitn, kappa, False, False, 1)
    else:
        if var_K:
            x_fitn[0]= [*Minimum[0][0],0]
        if var_Rho:
            x_fitn[1]= [*Minimum[0][1], 0]
        if var_w:
            x_fitn[2] = [*Minimum[0][2], 0]
        fitn_wert_y =  get_Wirkung(t, r, N, Intgrenze, T, *x_fitn, kappa, False, False, 1)

    return x_fitn, fitn_wert_y

def x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste):
    if variant ==1:
        x_fitn = [K_Liste, Rho_Liste, w_Liste]
    elif variant ==2:
        x_fitn = [K_Liste, Rho_Liste, []]
    elif variant ==3:
        x_fitn = [[], Rho_Liste, w_Liste]
    elif variant ==4:
        x_fitn = [K_Liste, [], w_Liste]
    elif variant ==5:
        x_fitn = [K_Liste, [], []]
    elif variant ==6:
        x_fitn = [[], Rho_Liste, []]
    elif variant ==7:
        x_fitn = [[],[],w_Liste]
    elif variant ==8:
        x_fitn = [[],[],[]]

    return x_fitn

def Variierer(N, K_randomi, rho_randomi,for_rho, w_randomi, variant, x_fitn,
        var_rho = True, var_K = True, var_w = True):
    if var_K:
        print('N', N)
        randomi2 = (2*np.random.random_sample(N) - 1)*(K_End - K_Anf)/10
        print('K_randomi and randomi2', K_randomi, randomi2)
        K_randomi = np.absolute(K_randomi + randomi2)
        print('K_randomi', K_randomi)
        x_fitn[0]= list(K_randomi)
    if var_rho:
        randomi2 = np.random.random()/10
        rho_randomi = rho_randomi + randomi2
        x_fitn[1] = get_rho_values(N,rho_randomi,for_rho,Constant,liste2, SameValues  =False)
    if var_w:
        print('N = ', N)
        randomi2 = np.random.random_sample(N)*(w_End - w_Anf)/10
        print('olaaaa', randomi2, w_randomi)
        w_randomi = w_randomi + randomi2
        x_fitn[2] = list(w_randomi)

    return x_fitn

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



def Minimierer(N, first, liste2, Minimum, var_K, var_Rho, var_w, K_Liste,
        Rho_Liste, w_Liste):
    variant = which_variant(var_K, var_Rho, var_w)

    x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
    x_fitn, fitn_wert_y = Zeilensparer(liste2, Minimum, N, first, x_fitn, var_K, var_Rho,var_w, variant)
    x_y_min = [x_fitn, fitn_wert_y]
    Mittelgr = 2
    K_M = x_y_min[1]/Mittelgr
    for i in range(Mittelgr):
        x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
        x_fitn, fitn_wert_x = Zeilensparer(liste2, Minimum, N, first, x_fitn,
                var_K, var_Rho, var_w, variant)

        if fitn_wert_y < x_y_min[1]:
            x_y_min[1] =fitn_wert_x
            x_y_min[0] = x_fitn
        K_M += fitn_wert_y/Mittelgr


    for m,tt in enumerate(temp):
        K_randomi = np.random.random_sample(N)*(K_End - K_Anf)
        rho_randomi = np.random.random()
        w_randomi = np.random.random_sample(N)*(w_End - w_Anf)
        if N==1 or N ==2:
            for_rho= 1
        else:
            for_rho = np.random.randint(N -1)
        for _ in range(4):
            x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
            x_fitn = Variierer (N, K_randomi, rho_randomi, for_rho
             ,w_randomi, variant, x_fitn, var_Rho, var_K, var_w)

            print('Rho_values2 = ', x_fitn[1])
            fitn_wert_y = get_Wirkung(t, r, N, Intgrenze, T,
                    *x_fitn, kappa, False, False, 1)
            boltzi = boltzmann(fitn_wert_x, fitn_wert_y, tt, K_M)
            if fitn_wert_x > fitn_wert_y:
                fitn_wert_x=fitn_wert_y
                if x_y_min[1] > fitn_wert_y:
                    x_y_min[0]= x_fitn
                    x_y_min[1]= fitn_wert_y

            elif (fitn_wert_x < fitn_wert_y) and (random.random() <= boltzi) :
                fitn_wert_x = fitn_wert_y
            print('Minimum = ', x_y_min[0], x_y_min[1])
    return x_y_min

if __name__ == "__main__":
    '''
    Variablen für das kausale System werden nachfolgend deklariert
    '''

    var_K, var_Rho, var_w = configfunktion('SollVarr?')
    K_Anf, K_End, K_List = configfunktion('Impuls')
    w_Anf, w_End, w_List = configfunktion('Frequenz')
    Constant, kappa, Rho_List = configfunktion('Einschraenkungen')
    Anzahl_N, first = configfunktion('Systemgroesse')


    T = np.pi
    r = si.symarray('r', 1)
    t = si.symarray('t', 1)

#   K_Anf = 0
#   K_End = 25
#   K_Liste = list(np.linspace(K_Anf, K_End, 10))


#   kappa = 0.0001
#   kappa_Anzahl = 1


#   #Rho_Liste = [0.8,0.07]#Rho_data.get_rho_values(N,Factor, Constant, liste2, SameValues  = True)

    x_Anf = 0
    x_End = np.pi

    Intgrenze = [x_Anf, x_End]
    Wirk = []

#   w_Anf = 0
#   w_End = 10

    '''
    Ab hier werden die fuer die Minimierung benötigten Variablen definiert
    '''
    Iter = 25                                      #Die Anzahl der temperaturiterationen
    hilfsarray_fuer_temp = np.linspace(0.01,5,Iter)
    Amplitude = 0.5                                #Amplitude der Fluktuation
                                                   #der Temperatur
    freq = np.pi                                    #Die Frequenz mit der die temp variieren soll
    halb = 0.001                                      #Der Halbwertswert fuer die
    temp = temperatur(Iter, hilfsarray_fuer_temp, halb, freq, Amplitude)
    print('temp =', temp)
    '''
    Ich passe das Minimieren zuerst fuer den eindim. Fall an, danach
    passe ich es für mehrere dim. an.
    '''
#   Anzahl_N = 4
    Minimum  = []

#   first = 3
    for j in range(first,Anzahl_N+1):
        w_Liste = eval(w_List)
        print('w_Liste', type(w_Liste))
        rho = symbols('rho0:j')
        liste = subs_coeffs(j, Constant, rho)
        print('liste', liste, j, Constant)
        liste2 = np.array(liste)/liste[0]
        print('liste2', liste2)



        #rho_randomi = j*0.1
        #Rho_Liste =get_rho_values(j,rho_randomi, Constant, liste, SameValues  =False)
        #print('Rho_Liste', Rho_Liste)
        K_Liste = eval(K_List)
        Constant = 1
        Factor = 1
        print('liste2', liste2)
        if var_Rho == True:
            Rho_Liste = get_rho_values(j, Factor, 1, Constant, liste2, SameValues  = True)
            print('Rho', Rho_Liste)

        else:
            Rho_Liste = eval(Rho_List)

        print('Rho_Liste = ', Rho_Liste)
        Minimum = Minimierer(j, first, liste2, Minimum, var_K,
                var_Rho, var_w, K_Liste, Rho_Liste, w_Liste )
        gg = open('Minimum3.txt', 'a')
        gg.write('Minimum fuer N = %d'%(j) + str(Minimum)+'\n')
        gg.close()
