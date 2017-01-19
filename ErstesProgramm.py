from __future__ import division
from sympy.utilities.lambdify import lambdify
from sympy.abc import r,t
from sympy.physics.quantum import TensorProduct
from pyx import *
from pyx.graph import axis
import scipy.integrate as integrate
import scipy.special,scipy.misc
import sympy as sy
import numpy as np
import random
import Rho_data

sigma1 = sy.Matrix([[0, 1],[1, 0]])
sigma2 = sy.Matrix([[0, -1j],[1j, 0]])
sigma3 = sy.Matrix([[1, 0],[0, -1]])
ident =  sy.Matrix([[1,0],[0,1]])

'''
Needed parts for the Integrand
'''

def prefactor(n):
    return (int(scipy.misc.factorial(n + 2))/
            (8 * sy.pi**(3/2)*sy.gamma('%d/%d'%(3+2*n,2))))

def diracEigenvalues(n):
    return (2*n + 1)/2

def integralKernelPlus(n, r, theta, phi):
    n=n-1
    lala1 = sy.jacobi(n, 1/2, 3/2, r)
    lala2 = sy.jacobi(n, 3/2, 1/2, r)
    return prefactor(n)*(sy.cos(r/2)*lala1*ident -
                        1j*sy.sin(r/2)*lala2*sigma_r(theta, phi))

def integralKernelMinus(n, r, theta, phi):
    n=n-1
    lala1 = sy.jacobi(n, 1/2, 3/2, r)
    lala2 = sy.jacobi(n, 3/2, 1/2, r)
    return prefactor(n)*(sy.cos(r/2)*lala1*ident +
                        1j*sy.sin(r/2)*lala2*sigma_r(theta, phi))

def sigma_r(theta, phi):
    return sy.sin(theta)*sy.cos(phi)*sigma1 + sy.sin(theta)*sy.sin(phi)*sigma2+ sy.cos(theta)*sigma3

def preMatrixPlus(n,K_Liste):
    a = K_Liste[n-1]
    b = sy.root(1 + a**2, 2)
    return sy.Matrix([[1- b, a],[-a, 1 + b]])

def preMatrixMinus(n,K_Liste):
    a = K_Liste[n-1]
    b = sy.root(1 + a**2, 2)
    return sy.Matrix([[1- b, -a],[ a, 1 + b]])

def projector(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    mat = sy.Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1,N + 1):
        mat += Rho_Liste[n-1]*sy.exp(-1j*w_Liste[n-1]*t)*(TensorProduct(preMatrixPlus(n,
            K_Liste),
            integralKernelPlus(n, r, theta, phi))
    + TensorProduct(preMatrixMinus(n,K_Liste),integralKernelMinus(n, r, theta,
        phi)))
    return mat

def projectorAdj(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    mat =  sy.Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1, N+1):
        mat += Rho_Liste[n-1]*sy.exp(1j*w_Liste[n-1]*t)*(TensorProduct(preMatrixPlus(n,K_Liste),
            integralKernelMinus(n, r, theta, phi))
            + TensorProduct(preMatrixMinus(n,K_Liste),integralKernelPlus(n, r, theta, phi)))
    return mat

def closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    return projector(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)*projectorAdj(t, r, theta,
            phi, N, Rho_Liste, w_Liste, K_Liste)

def lagrangian_without_bound_constr(t, r, theta, phi,N, Rho_Liste, w_Liste, K_Liste):
    sub = closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)
    return sy.trace(sub*sub) - sy.S(0.25) * sy.S(sy.trace(sub)**2)
'''
Needed parts for the Integrand

'''
'''
Constraints
'''

def boundedness_constraint(t,r,theta, phi, N, Rho_Liste, w_Liste, K_Liste, kappa):
    sub = closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)
    print ('kappa=', kappa)
    return sy.S(kappa)* sy.trace(sub)**2

'''
Constraints
'''

'''
Integrand and Action
'''
def integrand(t, point):
    return (max(lagr(t, point).real,0) + bound(t,point))*np.sin(point)**2

def get_Integrand_values2():
    y_Matrix =np.zeros(L)
    for j,point in enumerate(x_values):
        y_Matrix[j] = integrand(t, point)
    return y_Matrix

def get_Wirkung_fuer_kappa(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste, kappa):

    lagr = lambdify((t,r),
            ((lagrangian_without_bound_constr(t,r,0,0,N, Rho_Liste,
                w_Liste, K_Liste))), "numpy")

    bound = lambdify((t,r), (boundedness_constraint(t,r,0,0,N,Rho_Liste,
        w_Liste, K_Liste,kappa)), "numpy")

    integrand = lambda t1, r1 : (max(lagr(t1, r1),0).real +
            bound(t1,r1).real)*np.sin(r1)**2

    Wirkung = integrate.dblquad(lambda t1, r1 : (max(lagr(t1, r1),0).real +
            bound(t1,r1).real)*(np.sin(r1)**2)*np.exp(-(t1-T)**2/T)
    ,0,2*T, lambda r1: Intgrenze[0], lambda r1 : Intgrenze[1])[0]

    return Wirkung



#   def Hauptroutine_Integrand(x_Anf, x_End, K_Anf, K_End, K_Anzahl, L, kappa):
#       y_Matrix, x_Achse, Kurve_Names = get_Integrand_values(x_Anf, x_End, K_Anf,
#               K_End, K_Anzahl, L, kappa)

#       diag_plot(list(x_Achse), y_Matrix, 'Radius = r', 'Lagr(r)$*sin(r)^2$', Kurve_Names,'N=1_Integrand_kappa=%0.2f'%kappa, "br")

#   def Hauptroutine_Wirkung(Intgrenze, K_Liste, Rho_Liste, w_Liste, kappa):
#       y_Matrix, Kurve_Names, x_Achse = get_Wirkung_fuer_kappas(Intgrenze,
#               K_Liste, Rho_Liste, w_Liste, kappa)


#       diag_plot(list(x_Achse), y_Matrix, 'K', 'S(K)',Kurve_Names,'N=1_S(K)_kappaschar_K_Anzahl%d'%K_Anzahl, "br")



if __name__ == "__main__":
    t, theta, phi = symbols('t theta phi')

    T = 1 #Lebensdauer des Universums, wird fuer die Schwartzfunktion benoetigt

    L = 100

    K_Anzahl=10
    K_Anf = 0
    K_End = 2.5

    kappa = 0.01
    kappa_Anzahl = 1

    w_Liste = [0, 0]
    N = 1
    Rho_Liste = Rho_data.get_rho_values(N)

    x_Anf = 0
    x_End = np.pi

    x_values= np.linspace(x_Anf,x_End,L)
    y_Matrix =np.zeros((K_Anzahl, L))

    Kurve_Names=[]

    Intgrenze = [x_Anf, x_End]
    Wirk = []
    K_Liste =  list(np.linspace(K_Anf,K_End,K_Anzahl))
    for i in range(K_Anzahl):
        K2_Liste = [K_Liste[i]]
        #lagr = lambdify((t,r), lagrangian_without_bound_constr(t,r,0,0,N))
        #bound = lambdify((t,r),simplify(boundedness_constraint(t,r,0,0,N,kappa)))
        Wirk.append(get_Wirkung_fuer_kappa(t, r,N,Intgrenze, T,  K2_Liste, Rho_Liste, w_Liste,
            kappa))
        #Kurve_Names.append('k1='+'%3.2f'%K_Liste[1])
        #y_liste  = get_Integrand_values2()
        #y_Matrix[i,:] = y_liste

    Kurve_Names  = [kappa]
    diag_plot2(K_Liste, Wirk, 'K', 'Wirkung',
            Kurve_Names,'N=1_Wirkung_kappa=%0.2f'%kappa, "tr")

    #Hauptroutine_Integrand(x_Anf, x_End, K_Anf, K_End, K_Anzahl, L, kappa)
    #Hauptroutine_Wirkung(K_Anf, K_End, K_Anzahl, kappa_Anzahl)

