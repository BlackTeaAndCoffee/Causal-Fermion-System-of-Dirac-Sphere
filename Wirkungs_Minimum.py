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

def get_rho_values(N, Factor, for_rho, Constant, liste, SameValues = True):
    print('llist', liste)
    mix_list = np.arange(N-1)
    np.random.shuffle(mix_list)
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
            print('rho_val', rho_values[0], liste)
            zahl = wert - rho_values[0]
            if zahl == 0:
                return rho_values
            else:
                wert = zahl

            rho_values[N-1] = zahl/(diracEigenvalues(N)**2-0.25)
        else:
            for j,listval  in enumerate(liste):
                rho_values[for_rho] = Factor*(wert/liste[for_rho])
                zahl = wert - rho_values[for_rho]

                if zahl ==0:
                    return rho_values
                elif j == for_rho:
                    continue
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

def Zeilensparer(liste2, Minimum, N, x_fitn, var_K, var_Rho,var_w, variant):
    print('x_fitn', x_fitn, liste)
    N = j
    if j ==1:
        print('her  i am')
        if var_K:
            x_fitn[0]= np.random.random_sample(N)*(K_End - K_Anf)
        if var_Rho:
            rho_randomi = np.random.random()
            x_fitn[1]= get_rho_values(N,rho_randomi,0, Constant, liste2, SameValues  =False)
        if var_w:
            x_fitn[2] = np.random.random_sample(N)*(w_End - w_Anf)
        print('hola', x_fitn)
        fitn_wert_y =  get_Wirkung(t, r, N, Intgrenze, T, *x_fitn, kappa, False, False, 1)
    else:
        print(Minimum)
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

def Variierer(K_randomi, rho_randomi,for_rho, w_randomi, variant, x_fitn,
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
        randomi2 = np.random.random_sample(N)*(w_End - w_Anf)/10
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



def Minimierer(j, liste2, Minimum, var_K, var_Rho, var_w, K_Liste,
        Rho_Liste, w_Liste):
    variant = which_variant(var_K, var_Rho, var_w)
    N  = j
    x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
    x_fitn, fitn_wert_y = Zeilensparer(liste2, Minimum, j, x_fitn, var_K, var_Rho,var_w, variant)
    print('baboo', x_fitn)
    x_y_min = [x_fitn, fitn_wert_y]
    print(x_y_min)
    Mittelgr = 2
    print('x_y_min', x_y_min)
    K_M = x_y_min[1]/Mittelgr
    for i in range(Mittelgr):
        print('i',liste2)
        x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
        x_fitn, fitn_wert_x = Zeilensparer(liste2, Minimum, j, x_fitn,
                var_K, var_Rho, var_w, variant)

        if fitn_wert_y < x_y_min[1]:
            x_y_min[1] =fitn_wert_x
            x_y_min[0] = x_fitn
        K_M += fitn_wert_y/Mittelgr


    for m,tt in enumerate(temp):
        K_randomi = np.random.random_sample(N)*(K_End - K_Anf)
        print('K_randomi',K_randomi)
        rho_randomi = np.random.random()
        w_randomi = np.random.random_sample(N)*(w_End - w_Anf)
        if N==1 or N ==2:
            for_rho= 1
        else:
            for_rho = np.random.randint(N -1)
        print('for_rho', for_rho)
        for j in range(10):
            x_fitn = x_fitn_func(variant, K_Liste, Rho_Liste, w_Liste)
            x_fitn = Variierer (K_randomi, rho_randomi, for_rho
             ,w_randomi, variant, x_fitn, var_Rho, var_K, var_w)

            fitn_wert_y = get_Wirkung(t, r, N, Intgrenze, T,
                    *x_fitn, kappa, False, False, 1)
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
    return x_y_min

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
    K_Liste = list(np.linspace(K_Anf, K_End, 10))


    kappa = 0.0001
    kappa_Anzahl = 1

    Constant = 1
    rho = symbols('rho0:N')
    liste = subs_coeffs(N, Constant, rho)

    liste2 = np.array(liste)/liste[0]
    print('liste2', liste2)
    Rho_Liste = [0.8,0.07]#Rho_data.get_rho_values(N,Factor, Constant, liste2, SameValues  = True)

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
    Anzahl_N = 2
    Minimum  = []

    for j in range(1,Anzahl_N+1):
        #rho_randomi = j*0.1
        #Rho_Liste =get_rho_values(j,rho_randomi, Constant, liste, SameValues  =False)
        #print('Rho_Liste', Rho_Liste)
        K_Liste = [i for i in range(j)]
        Constant = j*0.1
        Factor = 1
        print('liste2', liste2)
        Rho_Liste = get_rho_values(j, Factor, 1, Constant, liste2, SameValues  = True)

        Minimum = Minimierer(j, liste2, Minimum, False,
                False, True, K_Liste, Rho_Liste, w_Liste )

