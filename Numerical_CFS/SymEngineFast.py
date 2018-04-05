from symengine import I
from sympy.printing import ccode
from Numerical_CFS.PyxPlot3d import *
from Numerical_CFS.configfunktion import configfunktion
from Numerical_CFS import get_data
import sympy as sy
import symengine as si
import numpy as np
from scipy import integrate
import random
import time
import os
import scipy.misc
import ctypes
import subprocess

sigma1 = si.zeros(2)
sigma1[1,0] = 1
sigma1[0,1] = 1
sigma2 = si.zeros(2)
sigma2[1,0] = I
sigma2[0,1] = -I
sigma3 = si.zeros(2)
sigma3[0,0] = 1
sigma3[1,1] = -1
ident = si.zeros(2)
ident[0,0] = 1
ident[1,1] = 1
r = si.symarray('r', 1)
t = si.symarray('t', 1)



'''
For more information on these Functions, look in to the master
thesis of Nikki Kilbertus, which is included on the
github site. If you are up to it, you can read from page 1, but from a
numerical perspective Chapter 2.3 Numerical Recipe (page 24ff) is important.

The main function is called get_Action.

Codewise speaking this programm can be split into Get_Integrand and Get_Action.
In Get_Integrand as the name already spills it out, i'm constructing the
Integrand, which will be later on integrated in Get_Action to get the Action.

Get_Integrand
=============

    Here my idea was to use sympy to get the analytic form of the
    integrand. So this way ideally only in the numeric integration
    i would get round off errors. The integration is done numerically
    because analytically it takes to much time. I have to do this integration
    a lot of times, for various parameters. So for smaller integrands it
    could be feasable to do the integration analytically but the integrand will
    get bigger and bigger for N getting  bigger. We are talking about several
    pages. So in comparison doing the integration numerically seems more
    plausible. I am not really able to say, how long the analytic integration
    would take for an expression, so i have no exlicit remarks on that.

    Another option would have been, to do everything with numpy arrays.
    I have done this once, if i remember correctly. This wasn't fast
    at all. The problem is, that i still have to construct the integrand
    out of huge amount of smaller terms. I think python as the overhead (i'm
    not sure if thats the correct name for it) causes it to be that slow.
    This code is in FullyNumeric.py.

Get_Action
==========

    In Get_Action i integrate over the Integrand i get via Get_Integrand.
    Depending on which method i decided to to use, different methods come
    to use.
    There is integrating with nquad from scipy, or nquad ctypes, or integrating
    with c.

'''
class C_F_S:
    '''
    :param N:    Codewise, this parameter determines the size of the Paramter-List\
                of the weights Rho, impulses K and frequencies w. \
                Physically speaking this is the Number of Shells of the Causal-Fermin system.
    :type N:     Integer.
    :param Integration_bound:   It's a List of two Lists. The two List should contain the\
                                the boundary constraints of the Double-Integral. First List\
                                should contain the time boundaries and the second should\
                                contain the space boundaries.
    :type Integration_bound:    List of two List, each containing two floats.
    :param T:    The life of the Universe, needed for the schwartzfunktion. It's\
                also the upper boundary of the time integration. 
    :type T:    float.
    :param K_Liste: Impulse variables for which the action is\
                calculated.
    :type K_Liste: List of N floats.
    :param w_Liste: Frequenc variables for which the action is\
                calculated.
    :type w_Liste: List of N floats.
    :param Rho_Liste: Weihts of the differents shells, for which the\
                action gets calculated.
    :type Rho_Lsite: List of N floats.
    :param kappa: It's needed for the boundednes constraint.
    :type kappa: float.
    :param Schwartzfunktion: For integer frequencies the schwartzfunktion can\
                be omitted. If so the time integration only needs to be done\
                over one period.
    :type Schartzfunktion: boolean.
    :param Comp_String:  If true opening, writing, closing of a file gets\
                skipped and the whole procedure of compiling and running gets\
                done in on line of code.
    :type Comp_String: boolean.
    :param Integration_Type: Integration_Type of integration. 1 for C-types, 2 for C, 3 for\
                testing , 4  Scipy-quadpack . It's not only the integration\
                method. According to this decision also the integrand gets\
                constructed.
    :type Integration_Type: 1,2,3 or 4.
    ''' 
    def __init__(self, N,  T, System_Parameters,Integration_bound = [[0, np.pi],[0, 2*np.pi]], Schwartzfunktion = True, 
                 Comp_String = False, Integration_Type = 1, Test_Action=False):
        self.N = N
        self.Integration_bound = Integration_bound
        self.T = T
        self.K_Liste = System_Parameters[0]
        self.Rho_Liste = System_Parameters[1]
        self.w_Liste = System_Parameters[2]
        self.kappa = System_Parameters[3]
        self.Schwartzfunktion = Schwartzfunktion
        self.Comp_String = Comp_String
        self.Integration_Type = Integration_Type
        self.Test_Action = Test_Action        

    def TensorProduct(self, mat11, mat22):
        """I needed a Tensorproduct for two matrices with symbolic elements.
        Other "Tensorproducts or Direct-Products" do not what i need. So
        i wrote this Product myself.

        :param mat11: Mathematically speaking a 2 by 2 matrix.
        :param mat22: Mathematically speaking a 2 by 2 matrix.
        :type mat11: 2 by 2 numpy.array.
        :type mat22: 2 by 2 numpy.array.
        :return: 4 by 4 matrix, which is the Tensorproduct.
        :rtype: 4 by 4 numpy.array.

        """
        mat1 = np.array(mat11, dtype= object).reshape(2,2)
        mat2 = np.array(mat22, dtype = object).reshape(2,2)
        return np.bmat([[mat1[0,0]*mat2, mat1[0,1]*mat2] ,[mat1[1,0]*mat2,
            mat1[1,1]*mat2]]).reshape(4,4)


    def prefactor(self, n):
        return int(scipy.misc.factorial(n + 2))/(8 *
                si.pi**(3/2)*si.gamma('%d/%d'%(6+2*n,2)))

    def diracEigenvalues(self, n):
        return (2*n + 1)/2

    def integralKernelPlus(self, n):
        n=n-1
        lala11 = sy.jacobi(n, 1/2, 3/2, r[0])
        lala21 = sy.jacobi(n, 3/2, 1/2, r[0])
        return self.prefactor(n)*(si.cos(r[0]/2)*lala11*ident -
                            I*si.sin(r[0]/2)*lala21*sigma3)

    def integralKernelMinus(self, n):
        n=n-1
        lala1 = sy.jacobi(n, 1/2, 3/2, r[0])
        lala2 = sy.jacobi(n, 3/2, 1/2, r[0])

        return self.prefactor(n)*(si.cos(r[0]/2)*lala1*ident +
                            I*si.sin(r[0]/2)*lala2*sigma3)

    def preMatrixPlus(self, a):
        b = si.sqrt(1 + a**2)
        matrix = si.zeros(2)

        matrix[0,0]= 1-b
        matrix[0,1]= a
        matrix[1,0]= -a
        matrix[1,1]= 1+b
        return matrix

    def preMatrixMinus(self, a):
        b = si.sqrt(1 + a**2)
        matrix = si.zeros(2)

        matrix[0,0]= 1-b
        matrix[0,1]= -a
        matrix[1,0]= a
        matrix[1,1]= 1+b
        return matrix

    def projector(self):
        mat = np.zeros((4,4), dtype = object)
        for n in range(1, self.N + 1):
            Koef =self.Rho_Liste[n-1]*si.exp(-I*self.w_Liste[n-1]*t[0])
            Term11 = self.TensorProduct(self.preMatrixPlus(self.K_Liste[n-1]),self.integralKernelPlus(n))
            Term21 = self.TensorProduct(self.preMatrixMinus(self.K_Liste[n-1]),self.integralKernelMinus(n))
            mat += Koef*(Term11 + Term21)
        return mat

    def projectorAdj(self):
        mat1 = np.zeros((4,4),  dtype = object)

        for n in range(1, self.N + 1):
            Koeff = self.Rho_Liste[n-1]*si.exp(I*self.w_Liste[n-1]*t[0])
            Term12 = self.TensorProduct(self.preMatrixPlus(self.K_Liste[n-1]),self.integralKernelMinus(n))
            Term22 = self.TensorProduct(self.preMatrixMinus(self.K_Liste[n-1]),self.integralKernelPlus(n))
            mat1 += Koeff*(Term12 +Term22)
        return mat1

    def closedChain(self):
        return np.dot(self.projector(),self.projectorAdj())

    def lagrangian_without_bound_constr(self):
        sub1 = self.closedChain()
        return np.trace(np.dot(sub1,sub1)) - 0.25 * np.trace(sub1)*np.trace(sub1)

    '''
    Constraints
    '''

    def boundadness_constraint(self):
        sub = self.closedChain()
        return self.kappa* np.trace(sub)**2

    '''
    Integrand and Action
    '''
    def Integrand(self):
        args = r[0] , t[0] 
        exprs = [self.lagrangian_without_bound_constr()]
        lagr = si.LambdifyCSE(args, exprs, real = False)#, as_scipy = True)
        exprs2 = [self.boundadness_constraint()]
        bound = si.Lambdify(args,exprs2, real = False)#, as_scipy = True)

        if self.Schwartzfunktion:
            integrand = 1/(2*np.pi**2) + (max(lagr(t1, r1).real,0) +
                    bound(t1,r1).real)*np.sin(r1)**2*np.exp(-(t1)**2/self.T)
        else:
            integrand = 1/(2*np.pi**2) + (max(lagr([t1, r1]).real,0) +
                    bound([t1,r1]).real)*np.sin(r1)**2

        return integrand

    def get_Integrand_with_c(self):
        '''
        For all of this i used the ccode function of sympy to generate C code for the integrand.

        What's happening is, that a C-File with the function, which will get
        integrated over later on, gets created and compiled.

        The generated C-code has not the form i need it to be, so some symbols
        need to be modified. Like exp to cexp and so on.

        I have an option self.Comp_String, which when True causes the C-File to
        get compiled in one OS-Command instead of seperate commands. I tried
        this out, because i thought that the first option could be faster.
        But it was not, if i remember correctly. Also that's not the Bottleneck
        of my programm, so i didn't focus on that any longer.

        Since the integration algorithm is in another file and is except for
        the integration boundaries fixed, i decided just to do it by hand every
        time. Those boundaries are for now also fixed, so it's ok.

        '''

        a = ccode(self.lagrangian_without_bound_constr())
        b = ccode(self.boundadness_constraint())
        g = open(get_data('funcprint2.txt'), 'r')
        g1 = g.read()
        g.close()
        Whole_String = ''
        Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
        Prototypen = 'static float xsav;\n'+ 'static float(*nrfunc)(float,float);\n'
        Integranddef = "float f(float r, float t)"+ "{return"
        Begin_von_Func = " 1/(2*cpow(M_PI,2)) + (fmax(creall("
        Func1 = a.replace("exp","cexp").replace("pow","cpow").replace("r_0","r").replace("t_0","t")
        Func1_End = "),0)"

        Whole_String += Bibliotheken + Prototypen + Integranddef + Begin_von_Func +Func1+ Func1_End

        Func2_Anf = "+ creall("
        Func2 =  b.replace("exp", "cexp").replace("pow","cpow" ).replace("r_0","r").replace("t_0","t")

        if self.Schwartzfunktion:
            Func22 = '))*sin(1.0L*r)*sin(1.0L*r)'
            Func2_End = '*cexp(-(cpow(t,2)/'+"cpow(%2.0f,2)))"%(self.T)+';'+'}\n'
            Whole_String += Func2_Anf + Func2 + Func22+ Func2_End + g1
        else:
            Func22 = '))*sin(1.0L*r)*sin(1.0L*r);}\n'
            Whole_String += Func2_Anf + Func2 + Func22+g1


        if self.Comp_String:
            os.system('gcc -o testlib2'+' << EOF '+Whole_String+ 'EOF -lm')
        else:
            f = open('testlib2.c', 'w')
            f.write(Whole_String)
            f.close()
            os.system('gcc -o testlib2 testlib2.c -lm')

    def get_Integrand_with_ctypes(self):

        a = ccode(self.lagrangian_without_bound_constr())
        b = ccode(self.boundadness_constraint())


        Whole_String = ''

        Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
        Integranddef = "double f(int n, double args[n])"+ "{return"
        Begin_von_Func = " 1/(2*cpow(M_PI,2)) + (fmax(creall(" #I added, 1/2*pi^2 because 
                                                               #the Integral can be zero,
                                                               #and then the integration takes 
                                                               #very long due to error control.
                                                               #By adding the factor above we make 
                                                               #the integral one bigger, and then 
                                                               #at the and we just need to substract it
                                                               #again.
        Func1 = a.replace("exp","cexp").replace("r_0","args[0]").replace("pow",
                "cpow").replace("t_0","args[1]")
        Func1_End = "),0)"
        Whole_String+= Bibliotheken + Integranddef + Begin_von_Func + Func1+Func1_End

        Func2_Anf = "+ creall("
        Func2 =  b.replace("exp", "cexp").replace("r_0","args[0]").replace("pow",
                "cpow").replace("t_0","args[1]")
        g = open(get_data('funcprint.c'), 'r')
        g1 = g.read()

        if self.Schwartzfunktion:
            Func22 = '))*sin(1.0L*args[0])*sin(1.0L*args[0])'
            Func2_End = '*cexp(-(cpow(args[1],2)/'+"cpow(%2.0f,2)))"%(self.T)+';'+'}\n'
            Whole_String += Func2_Anf + Func2 + Func22+ Func2_End + g1
        else:
            Func22 = '))*sin(1.0L*args[0])*sin(1.0L*args[0]);}\n'
            Whole_String += Func2_Anf + Func2 + Func22 + g1
        if self.Comp_String:
            os.system('gcc -x c -shared -o testlib2.so -fPIC'+' << EOF '+Whole_String+ 'EOF')
        else:

            f = open('testlib2.c', 'w')
            f.write(Whole_String)
            f.close()

            os.system('gcc -x c -shared -o testlib2.so -fPIC testlib2.c')

        g.close()

    def get_Test_Integrandt(self):
        '''
        This function is for testing the different integration routines,
        which are implemented here.

        Basically one can choose a function (The Analytic Result should best be
        known. Than one can compare the results. I have a sinus function, and
        and additional peak.

        So we can test, which routine is more capable in detecting the peak.)
        '''
        Whole_String = ''

        Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
        Integranddef = "double f(int n, double args[n])"+ "{return"
        Whole_String+= Bibliotheken + Integranddef

        g = open(get_data('funcprint.c'), 'r')
        g1 = g.read()

        Func2_1 ='(-1/(cpow(M_PI,0.5)*%2.3f))*'%(self.T)
        Func2_2 = "cexp(-(cpow(args[0]-1,2)/cpow(%2.3f,2)))"%(self.T)

        Func2_3 = '+  sin(args[0])'+';'+'}\n'
    #    Intdef2 = "double gf(int n, double args[n])"+ "{return 0. ;}"
    #    Intdef3 = "double hf(int n, double args[n])"+ "{return 2. ;}"

        Whole_String += Func2_1 + Func2_2 + Func2_3 + g1


        f = open('testlib2.c', 'w')
        f.write(Whole_String)
        f.close()

        os.system('gcc -x c -shared -o testlib2.so -fPIC testlib2.c')

        g.close()
    def control_Action(self):

        '''
        This function serves to test the minimizer.
        This is a higher dimensional
        parabola. In first order the action is also a
        parabola so if the minimizer won't work here
        it will definitely fail for the action.
        '''
        s = 0
        w_Min = [o for o in range(N)]
        R_Min = [o for o in range(N)]
        K_Min = [o for o in range(N)]

        Min = [w_Min, R_Min, K_Min]
        for ii in range(N):
            if var_K:
                s += (self.K_Liste[ii] - Min[0][ii])**2
            if var_Rho:
                s += (self.Rho_Liste[ii] - Min[1][ii])**2
            if var_w:
                s+= (self.w_Liste[ii] - Min[2][ii])**2
        return s



    def get_Action(self):

        '''
        I wanted to test out, how much faster it would be to do the integration with C.
        I tried out a simple double quad integration method (see integration2.c),
        but this simple method, wasn't able to detect delta peak like structures in
        the integrand.

        (I will include what i mean exactly with delta peak like structures. Because
        the more thinner the delta peak gets, at some point every method will
        fail.)

        So i would need to write a more refined integration method.
        At some point i will do that. But at first i wanted to try out, ctypes
        (I used ctypes and the nquad integration method of scipy, and indeed its
        much faster than simply nquad in scipy and also it detects much better the
        delta peaks than my simple quad integration method in C.) cython and maybe julia.

        '''
        if self.Test_Action:
            return self.control_Action()
        else:
            if self.Integration_Type == 1:
                '''Integration with Cytpes'''
                self.get_Integrand_with_ctypes()
                aa = time.time()

                lib=ctypes.CDLL('./testlib2.so')
                lib.f.restype = ctypes.c_double
                lib.f.argtypes = (ctypes.c_int,ctypes.c_double)
                zup = integrate.nquad(lib.f,[self.Integration_bound[0],self.Integration_bound[1]],
                        opts=[{'epsabs' :10e-8, 'epsrel': 10e-8 },
                            {'epsabs': 10e-8, 'epsrel' : 10e-8} ] )
                print('(Action, abserr)=',zup[0] -1, zup[1])
                handle = lib._handle # obtain the SO handle

                ctypes.cdll.LoadLibrary('libdl.so').dlclose(handle)

                tt = time.time()
                print('Passed time during integration in sec:',tt-aa)
                return zup[0] -1 #the 1 is due to the term i added in the integrand. 
                                 #It just cancels it out. 
            elif self.Integration_Type ==2:
                '''Integration with C. Up until here i basically construct the
                integrand  and then C takes over.'''#Stimmt was nicht.Also irgendwas
                self.get_Integrand_with_c()
                tt = time.time()

                result = subprocess.run(['./testlib2'], stdout=subprocess.PIPE)
                integr_val = result.stdout.decode('utf-8') - 1 #the 1 is due to the term i added in the integrand. 
                                 #It just cancels it out. 
 
                print('result', result.stdout.decode('utf-8')-1)
                aa = time.time()
                print('time', aa-tt)
                return float(integr_val)
            elif  self.Integration_Type ==3:
                '''Test_Typ'''
                self.get_Test_Integrandt(self.T)
                aa = time.time()
                lib=ctypes.CDLL('./testlib2.so')
                lib.f.restype = ctypes.c_double
                lib.f.argtypes = (ctypes.c_int,ctypes.c_double)
                result = integrate.nquad(lib.f, [self.Integration_bound[0]])
                tt = time.time()
                print(aa-tt)
                return result
            elif self.Integration_Type == 4:
                '''Scipy quadpack'''

                integrand = self.Integrand()
                tt = time.time()
                Action = integrate.nquad(lambda t1, r1 : integrand(t1,r1),
                        self.Integration_bound, opts=[{'epsabs' :10e-10, 'epsrel': 10e-10 },{
                            'epsabs': 10e-10, 'epsrel' : 10e-10} ])

                aa = time.time()
                print('Time it took to integrate in sec:',aa - tt)
                print(Action -1)
                return Action[0] -1#the 1 is due to the term i added in the integrand. 
                                 #It just cancels it out. 
 
            print('done')

        def get_integrand_values(self):
            '''
            This method is here for plotting the integrand.
            '''
            if self.Integration_Type == 1:
                self.get_Integrand_with_ctypes()
                os.system('gcc -o yolo testlib2.c -lm')
                os.system('./yolo')


            elif self.Integration_Type ==2:
                self.get_Integrand_with_c()

                os.system('gcc -o yolo testlib2.c -lm')
                os.system('./yolo')

            else: 
                raise ValueError("Integration_Type must be 1 or 2 for this to work.")


def Two_Dim_Pic():
    '''
    For N =2 i plotted the Action with respect to K_1 = K_Liste[1] and K_2
    '''
    K_Anzahl=N

    K_An = 25

    K1_Liste = np.linspace(0,K_End,K_An)
    K2_Liste = np.linspace(0,K_End,K_An)



    Kurve_Names=[]

    Integration_bound = [x_Anf, x_End]
    Wirk = []

    w_Liste = eval(w_List)
    K_Liste = eval(K_List)

    for jj in range(20):
        rho_1 = jj*0.05
        Rho_Liste = [rho_1, (1-rho_1)/3]

        norma = get_Action()
        d = open('NumbFor3d.txt', 'w')

        for k1 in K1_Liste:
            for k2 in K2_Liste:
                K_Liste=[k1, k2]
                Wert = get_Action()
                print(Wert)
                string = "%f %f %f \n" %(k1,k2,Wert/norma)
                d.write(string)
        d.close()
        PDFName =   "WirkungN%dRho1_%f_Rho2_%f_VarK1_K2_%05d" %(N,Rho_Liste[0],
                Rho_Liste[1], jj)
        Plot("NumbFor3d.txt", PDFName)
        f = open(PDFName +'.pdf', 'a')
        g = open("NumbFor3d.txt", "r")
        os.system("cat Settings.cfg >> " +PDFName + ".pdf")

        f.write('K_1, K_2, Action\n')
        os.system("cat NumbFor3d.txt >> " +PDFName+".pdf" )

def MainProg():
    var_K, var_Rho, var_w = configfunktion('Vary_Parameters')
    K_Anf, K_End, K_List = configfunktion('Impuls')
    w_Anf, w_End, w_List = configfunktion('Frequenz')
    Constant, kappa, Rho_List = configfunktion('Constraints')
    Anzahl_N, first, LifeTime = configfunktion('System_sizes')

    LifeTime = 2*np.pi #Liftime of the universe, it's needed for the Schwartzfunction
    N = 1

    kappa_Anzahl = 1

    x_Anf = 0.
    x_End = np.pi

    Kurve_Names=[]

    Integration_bound = [[x_Anf, x_End], [0, 1]]
    Wirk = []
    w_Liste = [0]#eval(w_List)
    K_Liste = [0]
    Rho_Liste = [1] #np.array([ 0.23596702, (1- 0.23596702)/3])#,0.1/6,0.1/10 ])#(0.1 +0.03)/6 ])
    Sys_Params = [K_Liste, Rho_Liste, w_Liste, kappa] 
    CFS_Action = C_F_S(N, Integration_bound, LifeTime, Sys_Params, Schwartzfunktion = False,  
               Comp_String = False, Integration_Type = 1)
    
    Wirkun = CFS_Action.get_Action()
    print(Wirkun)

if __name__ == "__main__":
    MainProg()
