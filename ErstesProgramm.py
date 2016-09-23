from __future__ import division
from sympy import *
from sympy.utilities.lambdify import lambdify
from sympy.abc import r,x
from sympy.physics.quantum import TensorProduct
import scipy.special
import scipy.misc
import numpy as np
import quartic as qp
import DiagXY as db

def Prefactor(n):
    return int(scipy.misc.factorial(n + 2))/(8 * pi**(2)*int(scipy.misc.factorial(n)))

def DiracEigenvalue(n):
    return (2*n + 1)/2

def IntegralKernelPlus(n, r, theta, phi):
    lala = jacobi(n, 1/2, 3/2, r)
    return Prefactor(n)*(cos(r/2)*lala*ident -
                        1j*sin(r/2)*lala*sigma_r(theta, phi))

def IntegralKernelMinus(n, r, theta, phi):
    lala = jacobi(n , 1/2, 3/2, r)
    return Prefactor(n)*(cos(r/2)*lala*ident +
                        1j*sin(r/2)*lala*sigma_r(theta, phi))

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
    Ich muss hier bei IntegralKernelMin und Plus bei den n's aufpassen, weil hier nach Nikki eigtl n-1 stehen
    sollte, aber halt bei der Schleife faengt n bei 0 statt 1 an, deswgen habe ich es mal so geeloest.
    '''

    mat = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1,N + 1):
        mat += rho[n-1]*exp(-1j*w[n-1]*t)*(TensorProduct(PreMatrixPlus(n), IntegralKernelPlus(n, r, theta, phi))
    + TensorProduct(PreMatrixMinus(n),IntegralKernelMinus(n, r, theta, phi)))
    return mat
def ProjectorAdj(t, r, theta, phi, N):
    mat =  Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    for n in range(1, N+1):
        mat += rho[n-1]*exp(1j*w[n-1]*t)*(TensorProduct(PreMatrixPlus(n),IntegralKernelMinus(n, r, theta, phi))
                + TensorProduct(PreMatrixMinus(n),IntegralKernelPlus(n, r, theta, phi)))
    return mat

def ClosedChain(t, r, theta, phi, N):
    return Projector(t, r, theta, phi, N)*ProjectorAdj(t, r, theta, phi, N)

def Lagrangian_without_Max(t, r, theta, phi,N):
    sub = ClosedChain(t, r, theta, phi, N)
    return trace(sub*sub) - 0.25 * trace(sub)**2

def TraceConstraint(N):
    add = 0
    for n in range(1, N+1):
        add += rho[n-1]*(DiracEigenvalue(n)**2 - 0.25)
    return add



if __name__ == "__main__":
    t, theta, phi = symbols('t theta phi')
    sigma1 = Matrix([[0, 1],[1, 0]])
    sigma2 = Matrix([[0, -1j],[1j, 0]])
    sigma3 = Matrix([[1, 0],[0, -1]])
    ident =  Matrix([[1,0],[0,1]])
    k = [10]
    w = [0]
    N = 1
    rho = [N]
    #simplify(ClosedChain(t, r, 0, 0, N))

    a1, a2, a3, a4, a5  = symbols('a1 a2 a3 a4 a5')


    a = simplify(ClosedChain(t, r, 0, 0, N))

    print simplify(Lagrangian_without_Max(t,r,0,0,N))
    lagr = lambdify((t,r), simplify(Lagrangian_without_Max(t,r,0,0,N)))

    L = 100
    x_values= []
    y_values =[]
    for i in range(L):
        ter = i*np.pi/L
        x_values.append(ter)
        yolo = np.sin(ter)**2*lagr(t,ter)
        print(yolo)
        y_values.append(yolo)

    db.f_gegn_x(list(x_values), y_values, 'Radius', 'Lagr', 'N=1_LAgr(r)')
    #Nullstellen = qp.x1plus(*c)
    #print Nullstellen
    #print simplify(Nullstellen)
    #print roots(simplify(a.charpoly(x).as_expr()),x)

