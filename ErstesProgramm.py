from __future__ import division
from sympy import *
from sympy.utilities.lambdify import lambdify
from sympy.abc import r,x
from sympy.physics.quantum import TensorProduct
from pyx import *
from pyx.graph import axis
import scipy.special
import scipy.misc
import numpy as np
import quartic as qp
import DiagXY as db


def Diag_plot(x_values, y_matrix, X_Title, Y_Title, PDFTitle):



    c = graph.graphxy(width=10,height=10,
                      x = graph.axis.linear(min = min(x_values), max = max(x_values),
                          title= X_Title),
                      y = graph.axis.linear(min = np.amin(y_matrix), max = np.amax(y_matrix),
                          title= Y_Title))

    for i in range(np.shape(y_Matrix)[0]):
        y_values = y_matrix[i,:]
        c.plot(graph.data.values(x = x_values, y = y_values),
                [graph.style.line([color.rgb.blue])])

    c.writePDFfile(PDFTitle)

def Prefactor(n):
    return int(scipy.misc.factorial(n + 2))/(8 * pi**(3/2)*gamma('%d/%d'%(3+2*n,2)))

def DiracEigenvalue(n):
    return (2*n + 1)/2

def IntegralKernelPlus(n, r, theta, phi):
    n=n-1
    lala1 = jacobi(n, 1/2, 3/2, r)
    lala2 = jacobi(n, 3/2, 1/2, r)
    return Prefactor(n)*(cos(r/2)*lala1*ident -
                        1j*sin(r/2)*lala2*sigma_r(theta, phi))

def IntegralKernelMinus(n, r, theta, phi):
    n=n-1
    lala1 = jacobi(n, 1/2, 3/2, r)
    lala2 = jacobi(n, 3/2, 1/2, r)
    return Prefactor(n)*(cos(r/2)*lala1*ident +
                        1j*sin(r/2)*lala2*sigma_r(theta, phi))

def sigma_r(theta, phi):
    return sin(theta)*cos(phi)*sigma1 + sin(theta)*sin(phi)*sigma2 + cos(theta)*sigma3

def PreMatrixPlus(n):
    a = k[n-1]
    b = root(1 + a**2, 2)
    return Matrix([[1- b, a],[-a, 1 + b]])

def PreMatrixMinus(n):
    a = k[n-1]
    b = root(1 + a**2, 2)
    return Matrix([[1- b, -a],[ a, 1 + b]])

def Projector(t, r, theta, phi, N):
    '''
    Ich muss hier bei IntegralKernelMin und Plus bei den n's aufpassen, weil hier nach
    Nikki eigtl n-1 stehen sollte, aber halt bei der Schleife faengt n bei 0 statt 1 an,
    deswgen habe ich es mal so geeloest.
    '''

    mat = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1,N + 1):
        mat += rho[n-1]*exp(-1j*w[n-1]*t)*(TensorProduct(PreMatrixPlus(n),
            IntegralKernelPlus(n, r, theta, phi))
    + TensorProduct(PreMatrixMinus(n),IntegralKernelMinus(n, r, theta, phi)))
    return mat
def ProjectorAdj(t, r, theta, phi, N):
    mat =  Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1, N+1):
        mat += rho[n-1]*exp(1j*w[n-1]*t)*(TensorProduct(PreMatrixPlus(n),
            IntegralKernelMinus(n, r, theta, phi))
            + TensorProduct(PreMatrixMinus(n),IntegralKernelPlus(n, r, theta, phi)))
    return mat

def ClosedChain(t, r, theta, phi, N):
    return Projector(t, r, theta, phi, N)*ProjectorAdj(t, r, theta, phi, N)

def Lagrangian_without_Max(t, r, theta, phi,N):
    sub = ClosedChain(t, r, theta, phi, N)
    return trace(sub*sub) - 0.25 * trace(sub)**2

def Lagrangian_without_Max2(t, r, theta, phi,N, kappa):
    sub = ClosedChain(t, r, theta, phi, N)
    return trace(sub*sub) +(kappa - 0.25) * trace(sub)**2

def TraceConstraint(N):
    add = 0
    for n in range(1, N+1):
        add += rho[n-1]*(DiracEigenvalue(n)**2 - 0.25)
    return add

def Lagrange_Schar(x_values, y_matrix, kappa):
    for m,k_value in enumerate(k_list):
        k = [k_value]
        a = simplify(ClosedChain(t, r, 0, 0, N))
        #lagr = lambdify((t,r), simplify(Lagrangian_without_Max(t,r,0,0,N)))
        lagr = lambdify((t,r), simplify( Lagrangian_without_Max2( t, r, 0,0,N, kappa)))
        for j,point in enumerate(x_values):
            yolo = np.sin(point)**2*lagr(t,point)

            y_Matrix[m,j]=yolo
    Diag_plot(list(x_values), y_Matrix, 'Radius', 'Lagr', 'N=1_LAgr(r)_kappa=%0.2f'%kappa)


if __name__ == "__main__":
    t, theta, phi = symbols('t theta phi')
    sigma1 = Matrix([[0, 1],[1, 0]])
    sigma2 = Matrix([[0, -1j],[1j, 0]])
    sigma3 = Matrix([[1, 0],[0, -1]])
    ident =  Matrix([[1,0],[0,1]])

    rho1, k1 = symbols('rho1 k1')
    w = [0]
    N = 1
    rho = [1]

    a1, a2, a3, a4, a5  = symbols('a1 a2 a3 a4 a5')

    L = 100

    K_Anzahl=10
    k_list=np.linspace(1,5,K_Anzahl)
    kappa = 0.01

    x_values= np.linspace(0,np.pi,L)
    y_Matrix =np.zeros((K_Anzahl, L))
    y_values=[]
    for m,k_value in enumerate(k_list):
        k = [k_value]
        a = simplify(ClosedChain(t, r, 0, 0, N))
        #lagr = lambdify((t,r), simplify(Lagrangian_without_Max(t,r,0,0,N)))
        lagr = lambdify((t,r), simplify( Lagrangian_without_Max2( t, r, 0,0,N, kappa)))
        for j,point in enumerate(x_values):
            yolo = np.sin(point)**2*lagr(t,point)

            y_Matrix[m,j]=yolo
    Diag_plot(list(x_values), y_Matrix, 'Radius', 'Lagr', 'N=1_LAgr(r)_kappa=%0.2f'%kappa)
