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

def diracEigenvalues(nn):
     return (2*nn + 1)/2


'''
Constraints
'''
def traceConstraint(N_t, Constant, rho):
    #dann benutze rho[0] usw.
    addd = Constant
    for n_t in range(1, N_t):
        addd -= rho[n_t-1]*(diracEigenvalues(n_t)**2 - 0.25)

    return addd


'''
Constraints
'''

def subs_coeffs(NN, Constant, rho):
    subs  = traceConstraint(NN, Constant, rho)
    liste = []
    if NN == 1:
        return [1]
    else:
        for ii in range(NN - 1):
            liste.append(Poly(subs).coeff_monomial(rho[ii]))
    return liste

def listmaker(NN, Constant):
    listim = []
    for ni in range(1, NN+1):
        listim.append((diracEigenvalues(ni)**2 -0.25)/2)

    return listim

def get_rho_values(N_r, Factor, for_rho, Constant, liste_in, SameValues = True):
    mix_list = np.arange(N_r)
    np.random.shuffle(mix_list)
    ultimo = mix_list[-1]
    while for_rho == ultimo:
        np.random.shuffle(mix_list)
        ultimo = mix_list[-1]
        print('haha')
    if SameValues:
        s_m = 0
        for il in range(1, N_r+1):
            s_m += (diracEigenvalues(il)**2 - 0.25)
        rho_values = [Constant/s_m for il in range(N_r)]
    else:
        wert = Constant
        rho_values = np.zeros(N_r)
        if N_r ==1:
            return [2]
#       elif N_r ==2:
#           subb = Factor*wert
#           rho_values[0] = subb/liste_in[0]
#           wert = wert - subb
#           if wert == 0:
#               return rho_values

#           rho_values[N_r-1] = 2*wert/(diracEigenvalues(N_r)**2-0.25)
        else:
            rho_values[for_rho] = Factor*wert/liste_in[for_rho]
            wert = wert - rho_values[for_rho]*liste_in[for_rho]
            if wert == 0:
                return rho_values


            for jl,listval  in enumerate(mix_list):
                if listval == ultimo:
                    continue
                elif listval == for_rho:
                    continue
                else:
                    print('listval',listval)
                    al = random.random()
                    rho_values[listval] = al*wert/liste_in[listval]
                    wert = wert - rho_values[listval]*liste_in[listval]
                    if wert < 0:
                        print('wert <0 ')
                    if wert == 0:
                        return rho_values

            rho_values[ultimo] = wert/liste_in[ultimo]
    s_t = 0
    print('rhovals, liste', rho_values, liste_in)
    for ll in range(len(rho_values)):
        s_t+= liste_in[ll]*rho_values[ll]
    print('Summe',s_t)
    return rho_values

def Zeilensparer(liste3, Minimum2, N_Z,first, x_fitn2, var_K, var_Rho,var_w,
        variant, Ausnahme):

    if Ausnahme:
        x_fitn2[0]= [*Minimum2[0],0]
        x_fitn2[1]= [*Minimum2[1], 0]
        x_fitn2[2] = [*Minimum2[2], 0]

        fitn_wert_y2 =  get_Wirkung(t, r, N_Z, Intgrenze, T, *x_fitn2, kappa, False, False, 1)
    elif first ==N_Z:
        if var_K:
            x_fitn2[0]= np.random.random_sample(N_Z)*(K_End - K_Anf)
        if var_Rho:
            rho_randomi = np.random.random()
            x_fitn2[1]= get_rho_values(N_Z,rho_randomi,0, Constant, liste3, SameValues  =False)
        if var_w:
            x_fitn2[2] = np.random.random_sample(N_Z)*(w_End - w_Anf)
        fitn_wert_y2 =  get_Wirkung(t, r, N_Z, Intgrenze, T, *x_fitn2, kappa, False, False, 1)
    else:
        if var_K:
            x_fitn2[0]= [*Minimum2[0][0],0]
        if var_Rho:
            x_fitn2[1]= [*Minimum2[0][1], 0]
        if var_w:
            x_fitn2[2] = [*Minimum2[0][2], 0]
        fitn_wert_y2 =  get_Wirkung(t, r, N_Z, Intgrenze, T, *x_fitn2, kappa, False, False, 1)
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

def Variierer(N_V, K_randomi, rho_randomi,for_rho, w_randomi, variant, x_fitn,
        var_rho = True, var_K = True, var_w = True):
    if var_K:
        randomi2 = (2*np.random.random_sample(N_V) - 1)*(K_End - K_Anf)/10
        K_randomi = np.absolute(K_randomi + randomi2)
        x_fitn[0]= list(K_randomi)
    if var_rho:
        randomi2 = np.random.random()/10
        rho_randomi = rho_randomi + randomi2
        x_fitn[1] = get_rho_values(N_V,rho_randomi,for_rho,Constant,liste2, SameValues  =False)
    if var_w:
        randomi2 = np.random.random_sample(N_V)*(w_End - w_Anf)/10
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



def Minimierer(N_M, first, liste_M, Minimum, var_K, var_Rho, var_w, K_Liste,
        Rho_Liste, w_Liste, Ausnahme):
    variant = which_variant(var_K, var_Rho, var_w)

    x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
    x_fitn, fitn_wert_y = Zeilensparer(liste_M, Minimum, N_M, first, x_fitn,
            var_K, var_Rho,var_w, variant,Ausnahme )
    x_y_min = [x_fitn, fitn_wert_y]
    Mittelgr = 2
    K_M = x_y_min[1]/Mittelgr
    for _ in range(Mittelgr):
        x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
        x_fitn, fitn_wert_x = Zeilensparer(liste_M, Minimum, N_M, first, x_fitn,
                var_K, var_Rho, var_w, variant, False)

        if fitn_wert_y < x_y_min[1]:
            x_y_min[1] =fitn_wert_x
            x_y_min[0] = x_fitn
        K_M += fitn_wert_y/Mittelgr


    for m,tt in enumerate(temp):
        K_randomi = np.random.random_sample(N_M)*(K_End - K_Anf)
        w_randomi = np.random.random_sample(N_M)*(w_End - w_Anf)
        rho_randomi = np.random.random()
        for_rho = np.random.randint(N_M)
        for _ in range(4):
            x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
            x_fitn = Variierer (N_M, K_randomi, rho_randomi, for_rho
             ,w_randomi, variant, x_fitn, var_Rho, var_K, var_w)

            fitn_wert_y = get_Wirkung(t, r, N_M, Intgrenze, T,
                    *x_fitn, kappa, False, False, 1)
            boltzi = boltzmann(fitn_wert_x, fitn_wert_y, tt, K_M)
            if fitn_wert_x > fitn_wert_y:
                fitn_wert_x=fitn_wert_y
                if x_y_min[1] > fitn_wert_y:
                    x_y_min[0]= x_fitn
                    x_y_min[1]= fitn_wert_y

            elif (fitn_wert_x < fitn_wert_y) and (random.random() <= boltzi) :
                fitn_wert_x = fitn_wert_y
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
    Ausnahme = configfunktion('Test')


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
    Iter = 30                                      #Die Anzahl der temperaturiterationen
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
    Minimum  = [[9.7903003749135973, 1.0428185324876262], np.array([ 0.24596702,
        0.25134433]), [0, 1]]
#   first = 3
    for SN in range(first,Anzahl_N+1):
        w_Liste = eval(w_List)
        print('w_Liste', type(w_Liste))
        liste2 = listmaker(SN, Constant)

        #rho_randomi = j*0.1
        #Rho_Liste =get_rho_values(j,rho_randomi, Constant, liste, SameValues  =False)
        #print('Rho_Liste', Rho_Liste)
        K_Liste = eval(K_List)
        Constant = 1
        Factor = 1
        print('liste2', liste2)
        if var_Rho == True:
            Rho_Liste = get_rho_values(SN, Factor, 1, Constant, liste2,
                    SameValues  = false)
            print('Rho', Rho_Liste)

        else:
            Rho_Liste = eval(Rho_List)

        print('Rho_Liste = ', Rho_Liste)
        Minimum = Minimierer(SN, first, liste2, Minimum, var_K,
                var_Rho, var_w, K_Liste, Rho_Liste, w_Liste, Ausnahme)
        gg = open('Minimum4.txt', 'a')
        gg.write('Minimum fuer N = %d'%(SN) + str(Minimum)+'\n')
        gg.close()
