from sympy.printing import ccode
from scipy import integrate
from symengine import I
from sympy import *
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

def subs_coeffs(N, Constant, rho):
    subs  = traceConstraint(N, Constant, rho)
    liste = []
    for i in range(N - 1):
        liste.append(Poly(subs).coeff_monomial(rho[i]))
    return liste

def get_rho_values(N, Factor, Constant, liste, SameValues = True):

    if SameValues:
        s = 0
        for i in range(1, N+1):
            s += (diracEigenvalues(i)**2 - 0.25)
        rho_values = [Constant/s for i in range(N)]
    else:
        wert = Constant
        rho_values = np.zeros(N)
        if N ==1:
            return [2]
        elif N ==2:
            rho_values[0] = Factor*(wert/liste[0])
            zahl = wert - rho_values[0]
            if zahl == 0:
                return rho_values
            else:
                wert = zahl

            rho_values[N-1] = zahl/(diracEigenvalues(N)**2-0.25)
            print('wadaaa', rho_values)
        else:
            for j in range(N-1):
                if j ==0:
                    rho_values[0] = Factor*(wert/liste[0])
                    zahl = wert - rho_values[0]
                    if zahl ==0:
                        rho_values[0]
                else:
                    a = random.random()
                    rho_values[j] = a* (wert/liste[j])    # hier steht ein Minus,weil die
                                                           # Koeffizienten negativ sind,
                                                           # finden aber als positive
                                                           # Zahlen Verwendung.
                    zahl = wert - rho_values[j]
                    if zahl == 0:
                        return rho_values
                    else:
                        wert = zahl

            rho_values[N-1] = zahl/(diracEigenvalues(N)**2-0.25)
    return rho_values

def Zeilensparer(liste):
    K_randomi = np.random.random_sample(N)*(K_End - K_Anf)
    rho_randomi = np.random.random()
    Rho_Liste = get_rho_values(N,rho_randomi, Constant, liste, SameValues  =False)
    print("yala", Rho_Liste)
    w_randomi = np.random.random_sample(N)*(w_End - w_Anf)

    fitn_wert_y =  get_Wirkung(t, r, N, Intgrenze, T, list(K_randomi),
            Rho_Liste, list(w_randomi), kappa, False, False, 1)
    x_fitn = [K_randomi, Rho_Liste, list(w_randomi)]
    return x_fitn, fitn_wert_y

def Variierer(K_randomi, rho_randomi, w_randomi, var_rho = True, var_K = True, var_w = True):
    if var_K:
        randomi2 = (2*np.random.random_sample(N) - 1)*(K_End - K_Anf)/10
        K_randomi = K_randomi + randomi2
        K_Liste = list(abs(K_randomi))
    if var_rho:
        randomi2 = np.random.random()/10
        rho_randomi = rho_randomi + randomi2
        Rho_Liste = get_rho_values(N,rho_randomi, Constant,liste, SameValues  =False)
    if var_w:
        randomi2 = np.random.random_sample(N)*(w_End - w_Anf)/10
        w_randomi = w_randomi + randomi2
        w_Liste = list(w_randomi)
    fitn_wert_y = get_Wirkung(t, r, N, Intgrenze, T, K_Liste, list(Rho_Liste), w_Liste, kappa, False, False, 1)
    x_fitn = [K_Liste, Rho_Liste, w_Liste]
    return fitn_wert_y, x_fitn

def RandVAlforVar(var_rho = True, var_K = True, var_w = True):
    if var_K:
        K_randomi = np.random.random_sample(N)*(K_End - K_Anf)
    elif var_Rho:
        rho_randomi = np.random.random()
    elif var_w:
        w_randomi = np.random.random_sample(N)*(w_End - w_Anf)

    return K_randomi, rho_randomi, w_randomi

if __name__ == "__main__":
    '''
    Variablen für das kausale System werden nachfolgend deklariert
    '''
    N = 2
    T = np.pi
    r = si.symarray('r', 1)
    t = si.symarray('t', 1)

    K_Anzahl=N
    K_Anf = 0
    K_End = 25

    kappa = 0.0001
    kappa_Anzahl = 1

    Constant = 1
    rho = symbols('rho0:N')
    liste = subs_coeffs(N, Constant, rho)

    liste2 = np.array(liste)/liste[0]

    #Rho_Liste = [0.8,0.07]#Rho_data.get_rho_values(N,Factor, Constant, liste2, SameValues  = True)

    x_Anf = 0
    x_End = np.pi

    Intgrenze = [x_Anf, x_End]
    Wirk = []

    w_Anf = 0
    w_End = 10
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

    x_fitn,fitn_wert_y = Zeilensparer(liste2)

    x_y_min = [x_fitn, fitn_wert_y]
    Mittelgr = 2
    K_M = x_y_min[1]/Mittelgr
    for i in range(Mittelgr):
        x_fitn, fitn_wert_x = Zeilensparer(liste2)

        if fitn_wert_y < x_y_min[1]:
            x_y_min[1] =fitn_wert_x
            x_y_min[0] = x_fitn
        K_M += fitn_wert_y/Mittelgr


    for m,tt in enumerate(temp):
        K_randomi = np.random.random_sample(N)*(K_End - K_Anf)
        rho_randomi = np.random.random()
        w_randomi = np.random.random_sample(N)*(w_End - w_Anf)

        for j in range(10):
            fitn_wert_y, x_fitn =  Variierer(K_randomi, rho_randomi, w_randomi, var_rho = True, var_K = True, var_w = True)

            boltzi = boltzmann(fitn_wert_x, fitn_wert_y, tt, K_M)
            print('boltzi = ', boltzi)
            if fitn_wert_x > fitn_wert_y:
                fitn_wert_x=fitn_wert_y
                if x_y_min[1] > fitn_wert_y:
                    x_y_min[0]= x_fitn
                    x_y_min[1]= fitn_wert_y

            elif (fitn_wert_x < fitn_wert_y) and (random.random() <= boltzi) :
                fitn_wert_x = fitn_wert_y
            print('Minimum = ', x_y_min[0], x_y_min[1])

