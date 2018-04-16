from Numerical_CFS.configfunktion import writeconfig
from Numerical_CFS.SymEngineFast import *
import mpmath as mp 
import symengine as si
import sympy as sy
import numpy as np
import scipy.misc as sm
tt = si.symarray('t', 8)
r = si.symarray('r', 1)
T = 1                   #float, Lifetime of the universe, which is needed for the
                        #       Schwartzfunktion
N = 1                   #integer, Shell-Number

kappa = 0            #float, Needed for boundary constant
kappa_Anzahl = 1        #float, Needed for plotting a diagram for various
                        # kappas

x_Anf = 0
x_End = np.pi
Integration_bound = [[x_Anf, x_End], [0,T]]

ind = 0                 #Which Impuls should be varied, 0 for K1 and 1 for K2


w_List = [0.]
K_List = [1]
Rho_List = [1] # i have to set this here for the initialing of the class C_FS
System_Parameters = [K_List, Rho_List, w_List, kappa]
CFS_Action = C_F_S(N, T, System_Parameters, Integration_bound, Schwartzfunktion = False, 
            Comp_String = False, Integration_Type = 1)


def zweitens(alpha, beta, n1):
    func = ((1 - r[0])**(alpha + n1))*((1 + r[0])**(beta+n1))
    return si.diff(func, *[r[0]]*n1)

def JacPol(alpha, beta, n1):
    first = ((-1)**n1/((2**n1)*sy.factorial(n1)))*((1 - r[0])**(-alpha))*((1 + r[0])**(-beta))*zweitens(alpha, beta, n1)
    return first

def Action_For_Test():
    K  = si.symbols('K')
    term1 = 1/(6*si.pi**8)
    term2 = (2*K*(15 + 31*K**2 + 9*K**4 + 9*K**6) + 3*(1 + K**2)**3*(5 - 3*K**2)*sy.acos(1 -2/(1 + K**2)))/((1 + K**2)**2)    
    funci = si.lambdify(K, [term1*term2])
    return funci
class TestClass():
    
    def test_TensorProdukt(self):
        a = np.array([[tt[0], tt[1]],[tt[2], tt[3]]], dtype = object).reshape(2,2)
        b = np.array([[tt[4], tt[5]],[tt[6], tt[7]]], dtype = object).reshape(2,2)

        TenMat = CFS_Action.TensorProduct(a, b)
        Result = np.array([[[ tt[0]*tt[4]],[tt[0]*tt[5]],[tt[1]*tt[4]],[tt[1]*tt[5]]],
                           [[ tt[0]*tt[6]],[tt[0]*tt[7]],[tt[1]*tt[6]],[tt[1]*tt[7]]],
                           [[ tt[2]*tt[4]],[tt[2]*tt[5]],[tt[3]*tt[4]],[tt[3]*tt[5]]],
                           [[ tt[2]*tt[6]],[tt[2]*tt[7]],[tt[3]*tt[6]],[tt[3]*tt[7]]]]).reshape(4,4)
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
    
    def test_JacobiPolys(self):
        n = 10
        m = 10
        listn = np.zeros(n, dtype = object)
        vals = np.zeros((n,m), dtype = object)
        for i in range(n):
            #print( sy.simplify(JacPol(sy.Rational(1,2),sy.Rational(3,2), i)),'=?=', sy.jacobi(i, 1/2,3/2, r[0]))
            listn[i] = si.lambdify(r[0], [sy.simplify(JacPol(sy.Rational(1,2),sy.Rational(3,2),i)) - 
                                       sy.jacobi(i, sy.Rational(1,2),sy.Rational(3,2), r[0])])   
            
        for i in range(n):
            for j in range(m):
                vals[i, j] = listn[i](j)[0]#, maxn = 10)

        print(vals) 
        np.testing.assert_almost_equal(vals, np.zeros((n,m)), 10)
 
    def test_integralKernelPlus(self):
        '''
        The idea in this test is to use the fact, that i know that if i 
        add the diagonal terms, the imaginary parts cancel each other out and
        the real terms add up. That's the idea for this and the following 3 tests. 
        Also this test is the most updated version. Here i do not use lambdify, 
        because of simplify, the terms really cancel to zero. Befor that i didn't 
        use simplify, so i had to use lambdify. I'm keeping the other tests the 
        way they are, so that i just can check how i did this. There were some 
        problems, i had to figure out and in future i maybe need to look it up.:) 
        '''
        n = 10
        m = 10
        listn = np.zeros(n, dtype = object)
        vals = np.zeros((n,m), dtype = object)
        for i in range(n):
            Plus0 = np.array(CFS_Action.integralKernelPlus(i+1))#, dtype = object).reshape(2,2)    
            print(Plus0)
            listn[i] =sy.simplify(-2*si.cos(r[0]/2)*CFS_Action.prefactor(i)*
                        JacPol(si.Rational(1,2),si.Rational(3,2), i) + Plus0[0,0] + 
                           Plus0[1,1])
        for i in range(n):
             vals[i] = listn[i]
        
        print(vals) 
        np.testing.assert_array_equal(listn, np.zeros(n), 10)
    
    def test_integralKernelPlus2(self):
        n = 10
        m = 10
        listn = np.zeros(n, dtype = object)
        vals = np.zeros((n,m), dtype = object)
        for i in range(n):
            Plus0 = np.array(CFS_Action.integralKernelPlus(i +1))#, dtype = object).reshape(2,2)    
            listn[i] =si.lambdify(r[0], [sy.simplify(Plus0[0,1] + Plus0[1,0])], real = False)
            
        for i in range(n):
            for j in range(m):
                vals[i, j] = listn[i](j)[0]
        
        print(vals) 
        np.testing.assert_almost_equal(vals, np.zeros((n,m)), 10)
    
    def test_integralKernelMinus(self):
        n = 10
        m = 10
        listn = np.zeros(n, dtype = object)
        vals = np.zeros((n,m), dtype = object)
        for i in range(n):
            Plus0 = np.array(CFS_Action.integralKernelMinus(i +1)).reshape(2,2)    
            print(Plus0)
            #Plus1 = np.array(CFS_Action.integralKernelMinus(i +1), dtype = object).reshape(2,2)
            listn[i] =si.lambdify(r[0], [sy.simplify(-2*si.cos(r[0]/2)*CFS_Action.prefactor(i)*
                              JacPol(si.Rational(1,2),si.Rational(3,2), i) + Plus0[0,0] + Plus0[1,1])], real = False)
            #print( sy.simplify(JacPol(1/2,3/2, i)),'=?=', sy.jacobi(i, 1/2,3/2, r[0]))
            #listn[i] =si.lambdify(r[0], [sy.simplify(JacPol(1/2,3/2, i)) - sy.jacobi(i, 1/2,3/2, r[0])])   
            
        for i in range(n):
            for j in range(m):
                vals[i, j] = listn[i](j)
        
        print(vals) 
        np.testing.assert_almost_equal(vals, np.zeros((n,m)), 10)
    
    def test_integralKernelMinus2(self):
        n = 10
        m = 10
        listn = np.zeros(n, dtype = object)
        vals = np.zeros((n,m), dtype = object)
        for i in range(n):
            Plus0 = np.array(CFS_Action.integralKernelMinus(i +1))#, dtype = object).reshape(2,2)    
            listn[i] =si.lambdify(r[0], [sy.simplify(Plus0[0,1] + Plus0[1,0])], real = False)
            
        for i in range(n):
            for j in range(m):
                vals[i, j] = listn[i](j)
        
        print(vals) 
        np.testing.assert_almost_equal(vals, np.zeros((n,m)), 10)

    def test_preMatrixPlus1(self):
        Matrix = CFS_Action.preMatrixPlus(r[0])        
        assert sy.simplify(Matrix[1,0] + Matrix[0,1]) == 0
    
    def test_preMatrixPlus2(self):
        Matrix = CFS_Action.preMatrixPlus(r[0])        
        assert sy.simplify(Matrix[0,0] + Matrix[1,1]) == 2

    def test_preMatrixMinus1(self):
        Matrix = CFS_Action.preMatrixPlus(r[0])        
        assert sy.simplify(Matrix[1,0] + Matrix[0,1]) == 0

    def test_preMatrixMinus2(self):
        Matrix = CFS_Action.preMatrixPlus(r[0])        
        assert sy.simplify(Matrix[0,0] + Matrix[1,1]) == 2
    
    def test_Lagrangian(self):
        """
        This Test only works for N  = 1, Rho=1.
        """
        n = 10
 
        LagVals = np.zeros(n)
        for i in range(n):
            KK = si.Rational(i,1)
            CFS_Action.K_Liste[0] = KK
            Lag = CFS_Action.lagrangian_without_bound_constr()
            delta = 1 - KK**2 + (1 + KK**2)*sy.cos(r[0])
            Res = si.pi**(-8)*8*(1 + KK**2)*si.cos(r[0]/2)**2*delta
            LagVals[i] = sy.simplify(Res- Lag)
        np.testing.assert_array_equal(LagVals, np.zeros(n))
    
    def test_Action(self):
        """
        For this test i used the analytic result of the action for N = 1 and 
        Rho_1 = 1 at page 32 in Nikkis Master-Thesis. 
        """
        n = 20
        Action = Action_For_Test()
        Listn = np.zeros(n)
        for i in range(n):
            Listn[i] = Action(i/10)
        print(Listn)
        Test = np.array([0.0008277342004438248, 0.0008309955105414712, 0.000840194562103519, 
                      0.0008536589176576792, 0.0008689658481976503, 0.0008835306967854006, 
                      0.0008951221231393858, 0.0009022625399987272, 0.0009042350910075936, 
                      0.0009010571906092057, 0.0008931757689961961, 0.0008812958950642769, 
                      0.0008662152758828479, 0.000848700978669209, 0.000829473187787505, 
                      0.0008091024283588821, 0.0007880785363976461, 0.0007667879751987183, 
                      0.0007455241133152157, 0.0007245150453482341])

        np.testing.assert_almost_equal(Listn, Test)


