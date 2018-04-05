from Numerical_CFS.configfunktion import writeconfig
from Numerical_CFS.SymEngineFast import *
import symengine as si
import sympy as sy
import numpy as np
import scipy.misc as sm
tt = si.symarray('t', 8)
r = si.symbols('r')
T = 1                   #float, Lifetime of the universe, which is needed for the
                        #       Schwartzfunktion
N = 2                   #integer, Shell-Number

kappa = 0            #float, Needed for boundary constant
kappa_Anzahl = 1        #float, Needed for plotting a diagram for various
                        # kappas

x_Anf = 0
x_End = np.pi
Integration_bound = [[x_Anf, x_End], [0,T]]

ind = 0                 #Which Impuls should be varied, 0 for K1 and 1 for K2


w_List = [0.,1.]
K_List = [0.,0.]
Rho_List = [1,0] # i have to set this here for the initialing of the class C_FS
System_Parameters = [K_List, Rho_List, w_List, kappa]
CFS_Action = C_F_S(N, Integration_bound, T, System_Parameters, Schwartzfunktion = False, 
            Comp_String = False, Integration_Type = 1)


def zweitens(alpha, beta, n):
    func = ((1 - r)**(alpha +n))*((1 + r)**(beta+n))
    return si.diff(func, *[r]*n)

def JacPol(alpha, beta, n):
    first = ((-1)**n/((2**n)*sm.factorial(n)))*((1 - r)**(-alpha))*((1 + r)**(-beta))*zweitens(alpha, beta, n)
    return first


class TestClass():
    
    def test_TensorProdukt(self):
        a = np.array([[tt[0], tt[1]],[tt[2], tt[3]]], dtype = object ).reshape(2,2)
        b = np.array([[tt[4], tt[5]],[tt[6], tt[7]]]).reshape(2,2)

        TenMat = CFS_Action.TensorProduct(a, b)
        Result = np.array([[[ tt[0]*tt[4]],[tt[0]*tt[5]],[tt[1]*tt[4]],[tt[1]*tt[5]]],
                           [[ tt[0]*tt[6]],[tt[0]*tt[7]],[tt[1]*tt[6]],[tt[1]*tt[7]]],
                           [[ tt[2]*tt[4]],[tt[2]*tt[5]],[tt[3]*tt[4]],[tt[3]*tt[5]]],
                           [[ tt[2]*tt[6]],[tt[2]*tt[7]],[tt[3]*tt[6]],[tt[3]*tt[7]]]], dtype = object).reshape(4,4)
        np.testing.assert_array_equal(TenMat,Result)

    def test_prefactor(self):
        n = 10
        to_test = np.zeros(n, dtype = object) 
        listn = np.zeros(n , dtype = object)
        for i in range(n):
            listn[i] = (sm.factorial(i + 2)*sm.factorial(i)*4**i)/(8*si.pi**2*(i+1/2)*sm.factorial(2*i))
            to_test[i] = CFS_Action.prefactor(i)
        np.testing.assert_almost_equal(sy.simplify(to_test - listn), np.zeros(n, dtype = object), decimal = 10)

    def test_diracEigenvalues(self):
        n = 10
        listn = np.zeros(n)
        for i in range(n):
            listn[i] = CFS_Action.diracEigenvalues(i)
        np.testing.assert_array_equal( np.array([1/2, 3/2, 5/2, 7/2, 9/2, 11/2, 13/2, 15/2, 17/2, 19/2]), listn)
    
    
    def test_integralKernelPlus(self):
        n = 10
        listn = np.zeros(n)
        for i in range(n):
            listn[i] =-2*CFS_Action.prefactor(i)*JacPol(1/2,2/3, i) + CFS_Action.integralKernelPlus(i +1)[0,0] + CFS_Action.integralKernelPlus(i +1)[1,1]
        np.testing.assert_almost_equal(listn, np.zeros(n))
            
