from __future__ import division
from sympy import *
from sympy.utilities.lambdify import lambdify
import numpy as np
import random

def diracEigenvalues(n):
     return (2*n + 1)/2


'''
Constraints
'''
def traceConstraint(N, Constant, rho):
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

def get_rho_values(N_r, Factor, for_rho, Constant, liste_in, SameValues = True):
    mix_list = np.arange(N_r-1)
    np.random.shuffle(mix_list)
    ultimo = mix_list[-1]
    while for_rho == ultimo:
        np.random.shuffle(mix_list)
        ultimo = mix_list[-1]
    if SameValues:
        s_m = 0
        for i in range(1, N_r+1):
            s_m += (diracEigenvalues(i)**2 - 0.25)
        rho_values = [Constant/s_m for i in range(N_r)]
    else:
        wert = Constant
        rho_values = np.zeros(N_r)
        if N_r ==1:
            return [2]

        else:
            rho_values[for_rho] = Factor*wert/liste_in[for_rho]
            wert = wert - rho_values[for_rho]*liste_in[for_rho]
            if wert ==0:
                return rho_values

            for jl,listval  in enumerate(mix_list):
                if listval == ultimo:
                    continue
                elif listval == for_rho:
                    continue
                else:
                    a = random.random()
                    rho_values[listval] = a*wert/liste_in[listval]    # hier steht ein Minus,weil die
                                                            # Koeffizienten negativ sind,
                                                           # finden aber als positive
                                                           # Zahlen Verwendung.
                    wert = wert - rho_values[listval]*liste_in[listval]
                    if wert < 0:
                        print('Wrong')
                    if wert == 0:
                        return rho_values

            rho_values[ultimo] = wert/liste_in[ultimo]
    s_t = 0

    for ll in range(len(rho_values)):
        s_t+= liste_in[ll]*rho_values[ll]
    #print('Summe',s_t)
    return rho_values

def listmaker(NN, Constant):
    listim = []
    for ni in range(1, NN+1):
        listim.append((diracEigenvalues(ni)**2 -0.25)/2)

    return listim


if __name__ == "__main__":
    N = 7
    Constant = 1
    for_rho = 0
    Factor = 1
    liste2 = listmaker(N, Constant)
    print('liste2', liste2, N, Constant)




    rho_values  = get_rho_values(N, 0.0000001, for_rho, Constant,liste2, False)
    print('jops',rho_values)

    for i in range(10000):
        Factor = np.random.random()
        rho_values  = get_rho_values(N, Factor, for_rho, Constant,liste2, False)
