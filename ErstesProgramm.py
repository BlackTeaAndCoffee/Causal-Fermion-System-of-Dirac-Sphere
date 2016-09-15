from numpy import linalg as LA
import scipy.integrate as integrate
import scipy.special
import scipy.misc
import numpy as np

def Prefactor(n):
    return scipy.misc.factorial(n + 2)/(8 * np.pi**(3/2)*scipy.special.gamma(n + 3/2))

def DiracEigenvalue(n):
    return (2*n + 1)/2

def IntegralKernelPlus(n, r, theta, phi):
    jacobi = scipy.special.jacobi(n, 1/2, 3/2)
    return Prefactor(n)*np.kron(ident, np.cos(r/2)*jacobi(r)*ident -
                        1j*np.sin(r/2)*jacobi(r)*sigma_r(theta, phi))

def IntegralKernelMinus(n, r, theta, phi):
    jacobi = scipy.special.jacobi(n, 1/2, 3/2)
    return Prefactor(n)*np.kron(ident, np.cos(r/2)*jacobi(r)*ident +
                        1j*np.sin(r/2)*jacobi(r)*sigma_r(theta, phi))

def sigma_r(theta, phi):
    return np.sin(theta)*np.cos(phi)*sigma1 + np.sin(theta)*np.sin(phi)*sigma2 + np.cos(theta)*sigma3

def PreMatrixPlus(n):
    a = k(n)
    b = np.sqrt(1 + a**2)
    return np.kron(np.array([[1- b, a],[-a, 1 + b]]), ident)

def PreMatrixMinus(n):
    a = k(n)
    b = np.sqrt(1 + a**2)
    return np.kron(np.array([[1- b, a],[ a, 1 + b]]), ident)

def Projector(t, r, theta, phi):
    mat = 0
    for n in range(N):
        mat += rho(n)*np.exp(-1j*w(n)*t)*(PreMatrixPlus(n)*IntegralKernelMinus(n-1, r, theta, phi)
                + PreMatrixMinus(n)*IntegralKernelPlus(n-1, r, theta, phi))
    return mat
def ProjectorAdj(t, r, theta, phi):
    mat = 0
    for n in range(N):
        mat += rho(n)*np.exp(1j*w(n)*t)*(PreMatrixPlus(n)*IntegralKernelMinus(n-1, r, theta, phi)
                + PreMatrixMinus(n)*IntegralKernelPlus(n-1, r, theta, phi))
    return mat

def ClosedChain(t, r, theta, phi):
    return np.dot(Projector(t, r, theta, phi), ProjectorAdj(t, r, theta, phi))

def Lagrangian(t, r, theta, phi):
    sub = ClosedChain(t, r, theta, phi)
    np.max(0 ,np.trace(np.dot(sub,sub )) - 0.25 * np.trace(sub)**2)

def TraceConstraint(N):
    add = 0
    for n in range(N):
        add += rho(n)*(DiracEigenvalue(n)**2 - 0.25)
    return add



if __name__ == "__main__":
    sigma1 = np.array([[0, 1],[1, 0]])
    sigma2 = np.array([[0, -1j],[1j, 0]])
    sigma3 = np.array([[1, 0],[0, -1]])
    ident = np.array([[1,0],[0,1]])
    k = [10]
    w = [0]
    N = 1
    rho = [1]
    EW, EV = LA.eig(ClosedChain())

