from sympy.printing import ccode
from scipy import integrate
from symengine import I
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


if __name__ == "__main__":
    '''
    Variablen für das kausale System werden nachfolgend deklariert
    '''
    N = 1
    T = np.pi
    r = si.symarray('r', 1)
    t = si.symarray('t', 1)

    K_Anzahl=N
    K_Anf = 0
    K_End = 2.5

    kappa = 0.003
    kappa_Anzahl = 1

    Rho_Liste = Rho_data.get_rho_values(N, SameValues  = True)

    x_Anf = 0
    x_End = np.pi

    Intgrenze = [x_Anf, x_End]
    Wirk = []


    w_Liste = [i for i in range(N)]

    '''
    Ab hier werden die fuer die Minimierung benötigten Variablen definiert
    '''
    Iter = 10                                      #Die Anzahl der temperaturiterationen
    hilfsarray_fuer_temp = np.linspace(0.01,50,Iter)
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
    x_fitn = random.random()*(K_End - K_Anf) #Startwert... random
    K_Liste = [x_fitn]
    fitn_wert_x = get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,kappa, False,False, 4)
    x_y_min = [x_fitn, fitn_wert_x]
    Mittelgr = 5
    K_M = x_y_min[1]/Mittelgr
    for i in range(Mittelgr):
        x_fitn = random.random()*(K_End - K_Anf) #Startwert... random
        print(x_fitn)
        K_Liste = [x_fitn]
        fitn_wert_y = get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,kappa, False,False, 4)
        if fitn_wert_y < x_y_min[1]:
            x_y_min[1] =fitn_wert_y
            x_y_min[0] = x_fitn
        K_M += fitn_wert_y/Mittelgr

    for m,tt in enumerate(temp):
        randomi = random.random()*(K_End - K_Anf)
        for j in range(10):
            randomi2 = (2*random.random() - 1)*(K_End - K_Anf)/10
            randomi = randomi + randomi2
            if K_Anf < randomi < K_End:

                K_Liste = [randomi]
                fitn_wert_y = get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste, kappa, False, False, 4)
                boltzi = boltzmann(fitn_wert_x, fitn_wert_y, tt, K_M)
                print('boltzi = ', boltzi)
                if fitn_wert_x > fitn_wert_y:
                    fitn_wert_x=fitn_wert_y
                    if x_y_min[1] > fitn_wert_y:
                        x_y_min[0]= randomi
                        x_y_min[1]= fitn_wert_y

                elif (fitn_wert_x < fitn_wert_y) and (random.random() <= boltzi) :
                    fitn_wert_x = fitn_wert_y
            print('Minimum = ', x_y_min[0], x_y_min[1])

