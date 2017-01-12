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
        print len(y_values)
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

'''
Needed parts for the Integrand
'''

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

def preMatrixPlus(n,K_Liste):
    a = K_Liste[n]
    b = root(1 + a**2, 2)
    return Matrix([[1- b, a],[-a, 1 + b]])

def preMatrixMinus(n,K_Liste):
    a = K_Liste[n]
    b = root(1 + a**2, 2)
    return Matrix([[1- b, -a],[ a, 1 + b]])

def projector(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    mat = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1,N + 1):
        mat += Rho_Liste[n-1]*exp(-1j*w_Liste[n-1]*t)*(TensorProduct(preMatrixPlus(n,
            K_Liste),
            integralKernelPlus(n, r, theta, phi))
    + TensorProduct(preMatrixMinus(n,K_Liste),integralKernelMinus(n, r, theta,
        phi)))
        print Rho_Liste, w_Liste, K_Liste, n
    return mat

def projectorAdj(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    mat =  Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1, N+1):
        mat += Rho_Liste[n-1]*exp(1j*w_Liste[n-1]*t)*(TensorProduct(preMatrixPlus(n,K_Liste),
            integralKernelMinus(n, r, theta, phi))
            + TensorProduct(preMatrixMinus(n,K_Liste),integralKernelPlus(n, r, theta, phi)))
    return mat

def closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    return projector(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)*projectorAdj(t, r, theta,
            phi, N, Rho_Liste, w_Liste, K_Liste)

def lagrangian_without_bound_constr(t, r, theta, phi,N, Rho_Liste, w_Liste, K_Liste):
    sub = closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)
    return trace(sub*sub) - 0.25 * trace(sub)**2
'''
Needed parts for the Integrand

'''
'''
Constraints
'''

def boundedness_constraint(t,r,theta, phi, N, Rho_Liste, w_Liste, K_Liste, kappa):
    sub = closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)
    print 'kappa=', kappa
    return kappa* trace(sub)**2

def traceConstraint():
    #dann benutze Rho_Liste[0] usw.
    add = 0
    for n in range(1, N):
        add += Rho_Liste[n-1]*(diracEigenvalues(n)**2 - 0.25)

    return add/(diracEigenvalues(N)**2 - 0.25)
'''
Constraints
'''
'''
Function fpr Integration
'''
def integration(funktion, AnfW, EndW):
    result = integrate.quad(lambda (t,r): funktion(t,r), AnfW, EndW)
    return result[0]
'''
Function fpr Integration
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

def get_Wirkung_fuer_kappa(Intgrenze, K_Liste, Rho_Liste, w_Liste, kappa):

    lagr = lambdify((t,r),
            simplify(lagrangian_without_bound_constr(t,r,0,0,N, Rho_Liste,
                w_Liste, K_Liste )))
    bound = lambdify((t,r),simplify(boundedness_constraint(t,r,0,0,N,Rho_Liste, w_Liste,
        K_Liste,kappa)))
    integrand = lambda r : (max(lagr(t, r),0) + bound(t,r))*np.sin(r)**2

    Wirkung = integrate.quad(integrand, Intgrenze[0],Intgrenze[1])[0]

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
    sigma1 = Matrix([[0, 1],[1, 0]])
    sigma2 = Matrix([[0, -1j],[1j, 0]])
    sigma3 = Matrix([[1, 0],[0, -1]])
    ident =  Matrix([[1,0],[0,1]])

    L = 100

    K_Anzahl=51
    K_Anf = 0
    K_End = 6

    kappa = 0
    kappa_Anzahl = 1

    Rho_Anzahl = 10
    Rho_Liste_Anf = 0
    Rho_Liste_End = 1

    Rho_Liste1, k1, k2 = symbols('Rho_Liste1 k1 k2')
    Rho_Liste = symbols('Rho_Liste0:N')

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
    K_Liste =  list(np.linspace(0,10,K_Anzahl))
    for i in range(K_Anzahl):
        K2_Liste = [0,K_Liste[i]]
        #lagr = lambdify((t,r), lagrangian_without_bound_constr(t,r,0,0,N))
        #bound = lambdify((t,r),simplify(boundedness_constraint(t,r,0,0,N,kappa)))
        Wirk.append(get_Wirkung_fuer_kappa(Intgrenze, K2_Liste, Rho_Liste, w_Liste,
            kappa))
        #Kurve_Names.append('k1='+'%3.2f'%K_Liste[1])
        #y_liste  = get_Integrand_values2()
        #y_Matrix[i,:] = y_liste

    Kurve_Names  = [kappa]
    diag_plot2(K_Liste, Wirk, 'K', 'Wirkung',
            Kurve_Names,'N=1_Wirkung_kappa=%0.2f'%kappa, "tr")

    #Hauptroutine_Integrand(x_Anf, x_End, K_Anf, K_End, K_Anzahl, L, kappa)
    #Hauptroutine_Wirkung(K_Anf, K_End, K_Anzahl, kappa_Anzahl)

