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
    print(subs)
    liste = []
    for i in range(N - 1):
        liste.append(Poly(subs).coeff_monomial(rho[i]))
    return liste

def get_rho_values(N, Factor, Constant, SameValues = True):

    if SameValues:
        s = 0
        for i in range(1, N+1):
            s += (diracEigenvalues(i)**2 - 0.25)
        rho_values = [Constant/s for i in range(N)]
    else:
        rho = symbols('rho0:N')
        liste = subs_coeffs(N, Constant, rho)
        print(liste)
        liste2 = np.array(liste)/liste[0]
        print(liste2)
        wert = Constant
        rho_values = np.zeros(N)
        if N ==1:
            return [2]
        else:
            for j in range(N-1):
                a = Factor*random.random()
                rho_values[j] = -a* (wert/liste[j])    # hier steht ein Minus,weil die
                                                       # Koeffizienten negativ sind,
                                                       # finden aber als positive
                                                       # Zahlen Verwendung.
                                                 #Das gefaellt mir nicht,wert durch
                                                 #liste[j] ist blödsinn. Siehe
                                                 #Wirkungs_Minummum.py

                zahl = wert - rho_values[j]
                if zahl == 0:
                    break
                else:
                    wert = zahl

            rho_values[N-1] = zahl/(diracEigenvalues(N)**2-0.25)
    return rho_values


if __name__ == "__main__":
    N = 5
    Constant = 1
    Factor = 1
    rho_values  = get_rho_values(N, Factor, Constant, False)
    print(rho_values)
'''
Attention:Das Programm ist nicht aktuell. Die akteuelle Version ist in
Wirkungs_Minimum.py!!!!
'''

