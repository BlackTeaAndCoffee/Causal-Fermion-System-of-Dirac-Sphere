from sympy import *
from sympy.utilities.lambdify import lambdify
from sympy.abc import r
from sympy.physics.quantum import TensorProduct
from pyx import *
from pyx.graph import axis
from ErstesProgramm import *
import scipy.integrate as integrate
import scipy.special,scipy.misc
import numpy as np
import random
import Rho_data

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
    N = 3
    Anzahl_an_Ber = 1
    T = np.pi

    t, theta, phi = symbols('t theta phi')

    K_Anzahl=10
    K_Anf = 0
    K_End = 2.5

    kappa = 0.01
    kappa_Anzahl = 1

    Rho_Liste = Rho_data.get_rho_values(N, SameValues  = True)
    print('rho : ', Rho_Liste)

    x_Anf = 0
    x_End = np.pi

    Intgrenze = [x_Anf, x_End]
    Wirk = []


    w_Liste = [-i for i in range(N)]
    K_Liste = [i for i in range(1,N+1)]

    print('w : ', w_Liste)
    print('k : ', K_Liste)

    for i in range(Anzahl_an_Ber):
        Wirk.append(get_Wirkung_fuer_kappa(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
            kappa, False))
    print(Wirk)

#   K_Liste =  list(np.linspace(K_Anf,K_End,K_Anzahl))
#   for i in range(K_Anzahl):
#       w_Liste = np.random.random_sample(N)
#       w_Liste[0] = 0

#       K2_Liste = [K_Liste[i]]
#       Wirk.append(get_Wirkung_fuer_kappa(t, r,N,Intgrenze, T,  K2_Liste, Rho_Liste, w_Liste,
#           kappa))

#   Kurve_Names  = [kappa]
#   diag_plot2(K_Liste, Wirk, 'K', 'Wirkung',
#           Kurve_Names,'N=1_Wirkung_kappa=%0.2f'%kappa, "tr")


