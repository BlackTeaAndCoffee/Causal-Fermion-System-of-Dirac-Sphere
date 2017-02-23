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
def traceConstraint(N,rho):
    #dann benutze rho[0] usw.
    add = 1.
    for n in range(1, N):
        add -= rho[n-1]*(diracEigenvalues(n)**2 - 0.25)

    return add


'''
Constraints
'''

def subs_coeffs(N, rho):
    subs  = traceConstraint(N, rho)
    liste = []
    for i in range(N - 1):
        liste.append(Poly(subs).coeff_monomial(rho[i]))
    return liste

def get_rho_values(N, SameValues = True):

    if SameValues:
        s = 0
        for i in range(1, N+1):
            s += (diracEigenvalues(i)**2 - 0.25)
        rho_values = [1/s for i in range(N)]

    else:
        rho = symbols('rho0:N')
        liste = subs_coeffs(N, rho)
        wert = 1
        rho_values = np.zeros(N)
        if N ==1:
            return [2]
        else:
            for j in range(N-1):
                a = random.random()
                rho_values[j] = -a* (wert/liste[j])    # hier steht ein Minus,weil die
                                                       # Koeffizienten negativ sind,
                                                       # finden aber als positive
                                                       # Zahlen Verwendung.
                zahl = wert - rho_values[j]
                if zahl == 0:
                    break
                else:
                    wert = zahl

            rho_values[N-1] = zahl/(diracEigenvalues(N)**2-0.25)
    return rho_values


if __name__ == "__main__":
    N = 1
    rho_values  = get_rho_values(N, True)



