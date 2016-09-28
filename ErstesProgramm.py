from __future__ import division
from sympy import *
from sympy.utilities.lambdify import lambdify
from sympy.abc import r,x
from sympy.physics.quantum import TensorProduct
from pyx import *
from pyx.graph import axis
import scipy.integrate as integrate
import scipy.special,scipy.misc
import numpy as np


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
        print len(y_values)
        dd.append(graph.data.values(x = x_Achse, y = y_values ,title = Kurv_Names[i]))

    c.plot(dd,[graph.style.line([color.gradient.Rainbow])])

    c.writePDFfile(PDFTitle)

def prefactor(n):
    return (int(scipy.misc.factorial(n + 2))/
            (8 * pi**(3/2)*gamma('%d/%d'%(3+2*n,2))))

def diracEigenvalues(n):
    return (2*n + 1)/2

def integralKernelPlus(n, r, theta, phi):
    n=n-1
    lala1 = jacobi(n, 1/2, 3/2, r)
    lala2 = jacobi(n, 3/2, 1/2, r)
    return prefactor(n)*(cos(r/2)*lala1*ident -
                        1j*sin(r/2)*lala2*sigma_r(theta, phi))

def integralKernelMinus(n, r, theta, phi):
    n=n-1
    lala1 = jacobi(n, 1/2, 3/2, r)
    lala2 = jacobi(n, 3/2, 1/2, r)
    return prefactor(n)*(cos(r/2)*lala1*ident +
                        1j*sin(r/2)*lala2*sigma_r(theta, phi))

def sigma_r(theta, phi):
    return sin(theta)*cos(phi)*sigma1 + sin(theta)*sin(phi)*sigma2 + cos(theta)*sigma3

def preMatrixPlus(n,k_value):
    a = k_value
    b = root(1 + a**2, 2)
    return Matrix([[1- b, a],[-a, 1 + b]])

def preMatrixMinus(n,k_value):
    a = k_value
    b = root(1 + a**2, 2)
    return Matrix([[1- b, -a],[ a, 1 + b]])

def projector(t, r, theta, phi, N,k_value):
    mat = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1,N + 1):
        mat += rho[n-1]*exp(-1j*w[n-1]*t)*(TensorProduct(preMatrixPlus(n,k_value),
            integralKernelPlus(n, r, theta, phi))
    + TensorProduct(preMatrixMinus(n, k_value),integralKernelMinus(n, r, theta, phi)))
    return mat

def projectorAdj(t, r, theta, phi, N,k_value):
    mat =  Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1, N+1):
        mat += rho[n-1]*exp(1j*w[n-1]*t)*(TensorProduct(preMatrixPlus(n,k_value),
            integralKernelMinus(n, r, theta, phi))
            + TensorProduct(preMatrixMinus(n,k_value),integralKernelPlus(n, r, theta, phi)))
    return mat

def closedChain(t, r, theta, phi, N,k_value):
    return projector(t, r, theta, phi, N,k_value)*projectorAdj(t, r, theta, phi, N,k_value)

def lagrangian_without_bound_constr(t, r, theta, phi,N,k_value):
    sub = closedChain(t, r, theta, phi, N,k_value)
    return trace(sub*sub) - 0.25 * trace(sub)**2

def boundedness_constraint(t,r,theta, phi, N,k_value, kappa):
    sub = closedChain(t, r, theta, phi, N,k_value)
    return kappa* trace(sub)**2

def traceConstraint(N):
    add = 0
    for n in range(1, N+1):
        add += rho[n-1]*(diracEigenvalues(n)**2 - 0.25)
    return add

def integration(funktion, AnfW, EndW):
    result = integrate.quad(lambda (t,r): funktion(t,r), AnfW, EndW)
    return result[0]

def get_Integrand_values(x_Anf, x_End, k_Anf, k_End, k_Anzahl, L, kappa):
    x_values= np.linspace(x_Anf,x_End,L)
    y_Matrix =np.zeros((K_Anzahl, L))
    k_liste =np.linspace(k_Anf,k_End,K_Anzahl)
    Kurve_Names =[]
    for m,k_value in enumerate(k_liste):
        Kurve_Names.append('k1='+'%3.2f'%k_value)
        lagr = lambdify((t,r), lagrangian_without_bound_constr(t,r,0,0,N,k_value))
        bound = lambdify((t,r),simplify(boundedness_constraint(t,r,0,0,N,k_value,kappa)))
        for j,point in enumerate(x_values):
            y_Matrix[m,j] =(max(lagr(t, point).real,0) + bound(t,point))*np.sin(point)**2
    return y_Matrix, x_values, Kurve_Names

def get_Wirkung_fuer_kappas(x_Anf, x_End, K_Anzahl, kappa_Anzahl):
    kappa_list = np.linspace(0,0.01,kappa_Anzahl)
    Kurve_Names = []
    x_values= np.linspace(0,np.pi,L)
    y_Matrix =np.zeros((kappa_Anzahl, K_Anzahl))

    k_liste =np.linspace(0,2.5,K_Anzahl)
    for m,kappa in enumerate(kappa_list):
        Kurve_Names.append('kappa='+'%5.4f'%kappa)
        print Kurve_Names
        for j,k_value in enumerate(k_liste):
            lagr = lambdify((t,r), simplify(lagrangian_without_bound_constr(t,r,0,0,N,k_value)))
            bound = lambdify((t,r),simplify(boundedness_constraint(t,r,0,0,N,k_value,kappa)))
            integrand = lambda r : (max(lagr(t, r),0) + bound(t,r))*np.sin(r)**2

            y_Matrix[m,j]= integrate.quad(integrand, 0,np.pi)[0]
    return y_Matrix, Kurve_Names, k_liste

def Hauptroutine_Integrand(x_Anf, x_End, K_Anf, K_End, K_Anzahl, L, kappa):
    y_Matrix, x_Achse, Kurve_Names = get_Integrand_values(x_Anf, x_End, K_Anf,
            K_End, K_Anzahl, L, kappa)

    diag_plot(list(x_Achse), y_Matrix, 'Radius = r', 'Lagr(r)$*sin(r)^2$', Kurve_Names,'N=1_Integrand_kappa=%0.2f'%kappa, "br")

def Hauptroutine_Wirkung(K_Anf, K_End, K_Anzahl, kappa_Anzahl):
    y_Matrix, Kurve_Names, x_Achse = get_Wirkung_fuer_kappas(K_Anf, K_End, K_Anzahl, kappa_Anzahl)


    diag_plot(list(x_Achse), y_Matrix, 'K', 'S(K)',Kurve_Names,'N=1_S(K)_kappaschar', "br")


if __name__ == "__main__":
    t, theta, phi = symbols('t theta phi')
    sigma1 = Matrix([[0, 1],[1, 0]])
    sigma2 = Matrix([[0, -1j],[1j, 0]])
    sigma3 = Matrix([[1, 0],[0, -1]])
    ident =  Matrix([[1,0],[0,1]])

    L = 100
    K_Anzahl=10
    K_Anf = 0
    K_End = 2.5
    kappa = 0.1

    kappa_Anzahl = 10
    rho1, k1 = symbols('rho1 k1')
    w = [0]
    N = 1
    rho = [1]

    x_Anf = 0
    x_End = np.pi


    Hauptroutine_Integrand(x_Anf, x_End, K_Anf, K_End, K_Anzahl, L, kappa)
#    y_Matrix, x_Achse, Kurve_Names = get_Integrand_values(x_Anf, x_End, K_Anf,
#            K_End, K_Anzahl, L, kappa)

#    diag_plot(list(x_Achse), y_Matrix, 'Radius = r', 'Lagr(r)$*sin(r)^2$', Kurve_Names,'N=1_Integrand_kappa=%0.2f'%kappa, "br")


    Hauptroutine_Wirkung(K_Anf, K_End, K_Anzahl, kappa_Anzahl)
 #    y_Matrix, Kurve_Names, x_Achse = get_Wirkung_fuer_kappas(K_Anf, K_End, K_Anzahl, kappa_Anzahl)


    #diag_plot(list(x_Achse), y_Matrix, 'K', 'S(K)',Kurve_Names,'N=1_S(K)_kappaschar', "br")
