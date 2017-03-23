from __future__ import division
from sympy.utilities.lambdify import lambdify
from sympy.abc import r,t
from sympy.physics.quantum import TensorProduct
from pyx import *
from pyx.graph import axis
from sympy import symbols
from sympy.printing import ccode
import scipy.integrate as integrate
import scipy.special,scipy.misc
import sympy as sy
import numpy as np
import random
import Rho_data
import os
import DiagXY
import ctypes

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
            K_Liste),integralKernelPlus(n, r, theta, phi))
    + TensorProduct(preMatrixMinus(n,K_Liste),integralKernelMinus(n, r, theta,
        phi)))
    return mat

def projectorAdj(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    mat =  sy.Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1, N+1):
        mat += Rho_Liste[n-1]*sy.exp(1j*w_Liste[n-1]*t)*(TensorProduct(preMatrixPlus(n,K_Liste),integralKernelMinus(n, r, theta, phi))
            + TensorProduct(preMatrixMinus(n,K_Liste),integralKernelPlus(n, r, theta, phi)))
    return mat

def closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):

    return projector(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)*projectorAdj(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)

def lagrangian_without_bound_constr(t, r, theta, phi,N, Rho_Liste, w_Liste, K_Liste):
    sub = closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)
    #print(sub)
    return sy.trace(sub*sub) - 0.25* sy.trace(sub)**2
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
def Integrand(t, r, N, Rho_Liste, w_Liste, K_Liste, kappa, T, Schwartzfunktion
        = True):
    lagr = lambdify((t,r),
            ((lagrangian_without_bound_constr(t,r,0,0,N, Rho_Liste,
                w_Liste, K_Liste))), "numpy")
    bound = lambdify((t,r), (boundedness_constraint(t,r,0,0,N,Rho_Liste,
        w_Liste, K_Liste,kappa)), "numpy")


    if Schwartzfunktion:
        integrand = lambda t1, r1 : (max(lagr(t1, r1).real,0) +
                bound(t1,r1).real)*np.sin(r1)**2*np.exp(-(t1)**2/T)
    else:
        integrand = lambda t1, r1 : (max(lagr(t1, r1).real,0) +
                bound(t1,r1).real)*np.sin(r1)**2

    return integrand

def get_Integrand_values2():
    y_Matrix =np.zeros(L)
    for j,point in enumerate(x_values):
        y_Matrix[j] = integrand(t, point)
    return y_Matrix

def get_Wirkung_fuer_kappa(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion):

    integrand = Integrand(t, r,N, Rho_Liste, w_Liste, K_Liste, kappa, T,
            Schwartzfunktion)

    Wirkung = integrate.nquad(lambda r1,t1 : integrand(t1,r1), [[0, np.pi],[0,2]])
#   Wirkung = integrate.dblquad(lambda r1,t1 :
#           integrand(t1, r1)
#   ,0,np.pi, lambda t1: 0, lambda t1 : 2, epsabs =
#   1.49e-16, epsrel = 1.49e-16)

    return Wirkung

def get_Integrand_with_c(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String):
    a = ccode(lagrangian_without_bound_constr(t,r,0,0,N, Rho_Liste,
                w_Liste, K_Liste))
    b = ccode(boundedness_constraint(t,r,0,0,N,Rho_Liste,w_Liste, K_Liste,kappa))
    print('ola')
    g = open('funcprint2.txt', 'r')
    g1 = g.read()
    g.close()
    GesamtString = ''
    Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
    Prototypen = 'static float xsav;\n'+ 'static float(*nrfunc)(float,float);\n'
    Integranddef = "float f(float r, float t)"+ "{return"
    Begin_von_Func = "(fmax(creall("
    Func1 = a.replace("exp","cexp").replace("pow","cpow")
    Func1_End = "),0)"

    GesamtString += Bibliotheken + Prototypen + Integranddef + Begin_von_Func +Func1+ Func1_End

    Func2_Anf = "+ creall("
    Func2 =  b.replace("exp", "cexp").replace("pow","cpow" )

    if Schwartzfunktion:
        Func22 = '))*sin(1.0L*r)*sin(1.0L*r)'
        Func2_End = '*cexp(-(cpow(t,2)/'+"cpow(%2.0f,2)))"%(T)+';'+'}\n'
        GesamtString += Func2_Anf + Func2 + Func22+ Func2_End + g1
    else:
        Func22 = '))*sin(1.0L*r)*sin(1.0L*r);}\n'
        GesamtString += Func2_Anf + Func2 + Func22+g1


    if Comp_String:
        os.system('gcc -o testlib2'+' << EOF '+GesamtString+ 'EOF -lm')
    else:
        f = open('testlib2.c', 'w')
        f.write(GesamtString)
        f.close()
        os.system('gcc -o testlib2 testlib2.c -lm')


def get_Integrand_with_ctypes(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String):
    a = ccode(lagrangian_without_bound_constr(t,r,0,0,N, Rho_Liste,
                w_Liste, K_Liste))
    b = ccode(boundedness_constraint(t,r,0,0,N,Rho_Liste,
        w_Liste, K_Liste,kappa))


    Gesamtstring = ''

    Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
    Integranddef = "double f(int n, double args[n])"+ "{return"
    Begin_von_Func = "(fmax(creall("
    Func1 = a.replace("exp","cexp").replace("r","args[0]").replace("pow", "cpow").replace("t","args[1]")
    Func1_End = "),0)"
    Gesamtstring+= Bibliotheken + Integranddef + Begin_von_Func + Func1+Func1_End

    Func2_Anf = "+ creall("
    Func2 =  b.replace("exp", "cexp").replace("r","args[0]").replace("pow", "cpow").replace("t","args[1]")
    g = open('funcprint.txt', 'r')
    g1 = g.read()

    if Schwartzfunktion:
        Func22 = '))*sin(1.0L*args[0])*sin(1.0L*args[0])'
        Func2_End = '*cexp(-(cpow(args[1],2)/'+"cpow(%2.0f,2)))"%(T)+';'+'}\n'
        Gesamtstring += Func2_Anf + Func2 + Func22+ Func2_End + g1
    else:
        Func22 = '))*sin(1.0L*args[0])*sin(1.0L*args[0]);}\n'
        Gesamtstring += Func2_Anf + Func2 + Func22 + g1

    if Comp_String:
        os.system('gcc -x c -shared -o testlib2.so -fPIC'+' << EOF '+Gesamtstring+ 'EOF')
    else:

        f = open('testlib2.c', 'w')
        f.write(Gesamtstring)
        f.close()

        os.system('gcc -x c -shared -o testlib2.so -fPIC testlib2.c')

    g.close()

def get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion = True, Comp_String = False, With_Ctypes= True):

    if With_Ctypes:
        get_Integrand_with_ctypes(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
            kappa, Schwartzfunktion, Comp_String)

        lib=ctypes.CDLL('/home/mustafa/Regensburg/Reproduktion_Von_Nikkis_Ergebnissen/Progs_mit_Sympy/testlib2.so')
        lib.f.restype = ctypes.c_double
        lib.f.argtypes = (ctypes.c_int,ctypes.c_double)
        print(integrate.nquad(lib.f,[[0,np.pi],[0,2]]))
    else:
        get_Integrand_with_c(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String)
        os.system('./testlib2')
    print('done')
def get_integrand_values(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion = True, Comp_String = False, With_Ctypes = True):

    if With_Ctypes:
        get_Integrand_with_ctypes(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
            kappa, Schwartzfunktion, Comp_String)
    else:
        get_Integrand_with_c(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String)

    os.system('gcc -o yolo testlib2.c -lm')
    os.system('./yolo')

if __name__ == "__main__":
    t, theta, phi = symbols('t theta phi')

    T = 1 #Lebensdauer des Universums, wird fuer die Schwartzfunktion benoetigt

    L = 100

    N = 5

    K_Anzahl=N
    K_Anf = 1
    K_End = 6

    kappa = 0.1
    kappa_Anzahl = 1

    w_Liste = [i for i in range(N)]
    Rho_Liste = Rho_data.get_rho_values(N)

    x_Anf = 0
    x_End = np.pi

    Kurve_Names=[]

    Intgrenze = [x_Anf, x_End]
    Wirk = []
    K_Liste =  list(np.linspace(K_Anf,K_End,K_Anzahl))

#    integrand1 = Integrand(t, r, N, Rho_Liste, w_Liste, K_Liste, kappa, T, Schwartzfunktion
#        = True)

#   xvalues = []
#   yvalues =[]

#   for i in range(100):
#      el = i*np.pi/100
#      yvalues.append(integrand1(1,el))
#      print (el, integrand1(1,el))
#      xvalues.append(el)
    #DiagXY.f_gegn_x(list(xvalues), yvalues,'r', 'lagr', 'clagrsimp' )
#   get_integrand_values(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
#       kappa, Schwartzfunktion = True, Comp_String = False)


    get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,kappa, True,
            False, True)
#   print(get_Wirkung_fuer_kappa(t, r, N, Intgrenze, T, K_Liste, Rho_Liste,
#       w_Liste,kappa, True))
