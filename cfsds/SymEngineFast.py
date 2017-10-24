from sympy.printing import ccode
from scipy import integrate
from symengine import I
from .PyxPlot3d import *
from .configfunktion import configfunktion
from . import get_data
import sympy as sy
import symengine as si
import numpy as np
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


'''
For more information on these Functions, look in to the master
thesis of Nikki Kilbertus, which is included on the
github site. If you are up to it, you can read from page 1, but from a
numerical perspective Chapter 2.3 Numerical Recipe (page 24ff) is important.

The main function is called get_Action.

Codewise speaking this programm can be split into Get_Integrand and Get_Action.
In Get_Integrand as the name already spills it out, i'm constructing the
Integrand, which will be later on integrated in Get_Action to get the Action.

*Get_Integrand
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

#Get_Action
    In Get_Action i integrate over the Integrand i get via Get_Integrand.
    Depending on which method i decided to to use, different methods come
    to use.
    There is integrating with nquad from scipy, or nquad ctypes, or integrating
    with c.

'''

def TensorProduct(mat11, mat22):
    '''
    Inputs are two matrices, and the output is the tensoproduct of those two.

    I needed a Tensorproduct for two matrices with symbolic elements.
    Other Tensorproducts oer Direct Products do not what i need.
    If there is a nicer way,
    please inform me.

    '''
    mat1 = np.array(mat11, dtype= object).reshape(2,2)
    mat2 = np.array(mat22, dtype = object).reshape(2,2)
    return np.bmat([[mat1[0,0]*mat2, mat1[0,1]*mat2] ,[mat1[1,0]*mat2,
        mat1[1,1]*mat2]]).reshape(4,4)


def prefactor(n):
    return int(scipy.misc.factorial(n + 2))/(8 *
            si.pi**(3/2)*sy.gamma('%d/%d'%(3+2*n,2)))

def diracEigenvalues(n):
    return (2*n + 1)/2

def integralKernelPlus(n, r):
    n=n-1
    lala11 = sy.jacobi(n, 1/2, 3/2, r[0])
    lala21 = sy.jacobi(n, 3/2, 1/2, r[0])
    return prefactor(n)*(si.cos(r[0]/2)*lala11*ident -
                        I*si.sin(r[0]/2)*lala21*sigma3)

def integralKernelMinus(n, r):
    n=n-1
    lala1 = sy.jacobi(n, 1/2, 3/2, r[0])
    lala2 = sy.jacobi(n, 3/2, 1/2, r[0])

    return prefactor(n)*(si.cos(r[0]/2)*lala1*ident +
                        I*si.sin(r[0]/2)*lala2*sigma3)

'''
def sigma_r():
    aa = si.sin(theta)*si.cos(phi)*sigma1 + si.sin(theta)*si.sin(phi)*sigma2
    bb =  si.cos(theta)*sigma3

    return aa + bb

In IntegralKernelPlus and Minus, at the end of the function is sigma3 instead
of sigma_r, because theta and phi are zero.

'''

def preMatrixPlus(a):
    b = si.sqrt(1 + a**2)
    matrix = si.zeros(2)

    matrix[0,0]= 1-b
    matrix[0,1]= a
    matrix[1,0]= -a
    matrix[1,1]= 1+b
    return matrix

def preMatrixMinus(a):
    b = si.sqrt(1 + a**2)
    matrix = si.zeros(2)

    matrix[0,0]= 1-b
    matrix[0,1]= -a
    matrix[1,0]= a
    matrix[1,1]= 1+b
    return matrix

def projector(t, r, N, Rho_Liste2, w_Liste2, K_Liste2):
    mat = np.zeros((4,4), dtype = object)
    for n in range(1,N + 1):
        Koef =Rho_Liste2[n-1]*sy.exp(-I*w_Liste2[n-1]*t[0])
        Term11 = TensorProduct(preMatrixPlus(K_Liste2[n-1]),integralKernelPlus(n, r))
        Term21 = TensorProduct(preMatrixMinus(K_Liste2[n-1]),integralKernelMinus(n, r))
        mat += Koef*(Term11 + Term21)
    return mat

def projectorAdj(t, r,  N, Rho_Liste3, w_Liste3, K_Liste3):
    mat1 = np.zeros((4,4),  dtype = object)

    for n in range(1, N+1):
        Koeff = Rho_Liste3[n-1]*sy.exp(I*w_Liste3[n-1]*t[0])
        Term12 = TensorProduct(preMatrixPlus(K_Liste3[n-1]),integralKernelMinus(n, r))
        Term22 = TensorProduct(preMatrixMinus(K_Liste3[n-1]),integralKernelPlus(n, r))
        mat1 += Koeff*(Term12 +Term22)
    return mat1

def closedChain(t, r,  N, Rho_Liste, w_Liste, K_Liste):
    return np.dot(projector(t, r,  N, Rho_Liste, w_Liste,
        K_Liste),projectorAdj(t, r, N, Rho_Liste, w_Liste, K_Liste))

def lagrangian_without_bound_constr(t, r, N, Rho_Liste, w_Liste, K_Liste):
    sub1 = closedChain(t, r,  N, Rho_Liste, w_Liste, K_Liste)
    return np.trace(np.dot(sub1,sub1)) - 0.25 * np.trace(sub1)*np.trace(sub1)

'''
Constraints
'''

def boundedness_constraint(t,r, N, Rho_Liste, w_Liste, K_Liste, kappa):
    sub = closedChain(t, r,  N, Rho_Liste, w_Liste, K_Liste)
    return kappa* np.trace(sub)**2

'''
Integrand and Action
'''
def Integrand(t, r, N, Rho_Liste, w_Liste, K_Liste, kappa, T, Schwartzfunktion
        = True):
    args = np.concatenate((t,r))
    exprs = [lagrangian_without_bound_constr(t,r,N,Rho_Liste,w_Liste,K_Liste)]
    lagr = si.Lambdify(args, exprs, real = False)
    exprs2 = [boundedness_constraint(t,r,N,Rho_Liste,w_Liste, K_Liste,kappa)]
    bound = si.Lambdify(args,exprs2, real = False)

    if Schwartzfunktion:
        integrand = lambda t1, r1: (max(lagr(t1, r1).real,0) +
                bound(t1,r1).real)*np.sin(r1)**2*np.exp(-(t1)**2/T)
    else:
        integrand = lambda t1, r1 : (max(lagr([t1, r1]).real,0) +
                bound([t1,r1]).real)*np.sin(r1)**2

    return integrand

def get_Integrand_with_c(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String):
    '''
    For all of this i used the ccode function of sympy to generate C code for the integrand.

    What's happening is, that a C-File with the function, which will get
    integrated over later on, gets created and compiled.

    The generated C-code has not the form i need it to be, so some symbols
    need to be modified. Like exp to cexp and so on.

    I have an option Comp_String, which when True causes the C-File to
    get compiled in one OS-Command instead of seperate commands. I tried
    this out, because i thought that the first option could be faster.
    But it was not, if i remember correctly. Also that's not the Bottleneck
    of my programm, so i didn't focus on that any longer.

    Since the integration algorithm is in another file and is except for
    the integration boundaries fixed, i decided just to do it by hand every
    time. Those boundaries are for now also fixed, so it's ok.

    '''

    a = ccode(lagrangian_without_bound_constr(t,r,N, Rho_Liste,
                w_Liste, K_Liste))
    b = ccode(boundedness_constraint(t,r,N,Rho_Liste,w_Liste, K_Liste,kappa))
    g = open(get_data('funcprint2.txt'), 'r')
    g1 = g.read()
    g.close()
    Whole_String = ''
    Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
    Prototypen = 'static float xsav;\n'+ 'static float(*nrfunc)(float,float);\n'
    Integranddef = "float f(float r, float t)"+ "{return"
    Begin_von_Func = "(fmax(creall("
    Func1 = a.replace("exp","cexp").replace("pow","cpow").replace("r_0","r").replace("t_0","t")
    Func1_End = "),0)"

    Whole_String += Bibliotheken + Prototypen + Integranddef + Begin_von_Func +Func1+ Func1_End

    Func2_Anf = "+ creall("
    Func2 =  b.replace("exp", "cexp").replace("pow","cpow" ).replace("r_0","r").replace("t_0","t")

    if Schwartzfunktion:
        Func22 = '))*sin(1.0L*r)*sin(1.0L*r)'
        Func2_End = '*cexp(-(cpow(t,2)/'+"cpow(%2.0f,2)))"%(T)+';'+'}\n'
        Whole_String += Func2_Anf + Func2 + Func22+ Func2_End + g1
    else:
        Func22 = '))*sin(1.0L*r)*sin(1.0L*r);}\n'
        Whole_String += Func2_Anf + Func2 + Func22+g1


    if Comp_String:
        os.system('gcc -o testlib2'+' << EOF '+Whole_String+ 'EOF -lm')
    else:
        f = open('testlib2.c', 'w')
        f.write(Whole_String)
        f.close()
        os.system('gcc -o testlib2 testlib2.c -lm')

def get_Integrand_with_ctypes(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String):

    a = ccode(lagrangian_without_bound_constr(t,r,N, Rho_Liste,
                w_Liste, K_Liste))
    b = ccode(boundedness_constraint(t,r,N,Rho_Liste,
        w_Liste, K_Liste,kappa))


    Whole_String = ''

    Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
    Integranddef = "double f(int n, double args[n])"+ "{return"
    Begin_von_Func = "(fmax(creall("
    Func1 = a.replace("exp","cexp").replace("r_0","args[0]").replace("pow",
            "cpow").replace("t_0","args[1]")
    Func1_End = "),0)"
    Whole_String+= Bibliotheken + Integranddef + Begin_von_Func + Func1+Func1_End

    Func2_Anf = "+ creall("
    Func2 =  b.replace("exp", "cexp").replace("r_0","args[0]").replace("pow",
            "cpow").replace("t_0","args[1]")
    g = open(get_data('funcprint.c'), 'r')
    g1 = g.read()

    if Schwartzfunktion:
        Func22 = '))*sin(1.0L*args[0])*sin(1.0L*args[0])'
        Func2_End = '*cexp(-(cpow(args[1],2)/'+"cpow(%2.0f,2)))"%(T)+';'+'}\n'
        Whole_String += Func2_Anf + Func2 + Func22+ Func2_End + g1
    else:
        Func22 = '))*sin(1.0L*args[0])*sin(1.0L*args[0]);}\n'
        Whole_String += Func2_Anf + Func2 + Func22 + g1

    if Comp_String:
        os.system('gcc -x c -shared -o testlib2.so -fPIC'+' << EOF '+Whole_String+ 'EOF')
    else:

        f = open('testlib2.c', 'w')
        f.write(Whole_String)
        f.close()

        os.system('gcc -x c -shared -o testlib2.so -fPIC testlib2.c')

    g.close()

def get_Test_Integrandt(T):
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

    Func2_1 ='(-1/(cpow(M_PI,0.5)*%2.3f))*'%(T)
    Func2_2 = "cexp(-(cpow(args[0]-1,2)/cpow(%2.3f,2)))"%(T)

    Func2_3 = '+  sin(args[0])'+';'+'}\n'
#    Intdef2 = "double gf(int n, double args[n])"+ "{return 0. ;}"
#    Intdef3 = "double hf(int n, double args[n])"+ "{return 2. ;}"

    Whole_String += Func2_1 + Func2_2 + Func2_3 + g1


    f = open('testlib2.c', 'w')
    f.write(Whole_String)
    f.close()

    os.system('gcc -x c -shared -o testlib2.so -fPIC testlib2.c')

    g.close()


def get_Action(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion = True, Comp_String = False, Type = 1):

    '''
    Inputs are

    t,r          both 1 dimensional symarrays, these are also the integration
                 variables.

    N            Integer, from 1,2,....  .
                 Physically speaking the Shell-Number of the Causal-Fermin system

    Integration_bound   It's a float number

    T            Float, and the life of the Universe, needed for the schwartzfunktion

    K_Liste      List of floats, Impulse variables for which the action is
                 calculated

    w_Liste      List of floats, Frequenc variables for which the action is
                 calculated

    Rho_Liste    List of floats, Weihts of the differents shells, for which the
                 action gets calculated

    kappa        Float, it's needed for the boundednes constraint

    Schwartzfunktion boolean, For integer frequencies the schwartzfunktion can
                 be omitted, because then the time integration will not be
                 necessary

    Comp_String  boolean, if true opening, writing, closing of a file gets
                 skipped and the whole procedure of compiling and running gets
                 done in on line of code.

    Type         integer, type of integration. 1 for C-types, 2 for C, 3 for
                 testing , 4  Scipy-quadpack . It's not only the integration
                 method. According to this decision also the integrand gets
                 constructed.

    '''



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

    if Type == 1:
        '''Integration with Cytpes'''
        get_Integrand_with_ctypes(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
            kappa, Schwartzfunktion, Comp_String)
        aa = time.time()

        lib=ctypes.CDLL('./testlib2.so')
        lib.f.restype = ctypes.c_double
        lib.f.argtypes = (ctypes.c_int,ctypes.c_double)
        zup = integrate.nquad(lib.f,[Integration_bound[0],Integration_bound[1]],
                opts=[{'epsabs' :10e-8, 'epsrel': 10e-8 },
                    {'epsabs': 10e-10, 'epsrel' : 10e-10} ] )
        print('(Action, abserr)=',zup)
        handle = lib._handle # obtain the SO handle

        ctypes.cdll.LoadLibrary('libdl.so').dlclose(handle)

        tt = time.time()
        print('Passed time during integration in sec:',tt-aa)
        return zup[0]
    elif Type ==2:
        '''Integration with C. Up until here i basically construct the
        integrand  and then C takes over.'''#Stimmt was nicht.Also irgendwas
        get_Integrand_with_c(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String)
        tt = time.time()

        result = subprocess.run(['./testlib2'], stdout=subprocess.PIPE)
        integr_val = result.stdout.decode('utf-8')
        print('result', result.stdout.decode('utf-8'))
        aa = time.time()
        print('time', aa-tt)
        return float(integr_val)
    elif  Type ==3:
        '''Test_Typ'''
        get_Test_Integrandt(T)
        aa = time.time()
        lib=ctypes.CDLL('./testlib2.so')
        lib.f.restype = ctypes.c_double
        lib.f.argtypes = (ctypes.c_int,ctypes.c_double)
        result = integrate.nquad(lib.f, [[0,np.pi]])
        tt = time.time()
        print(aa-tt)
        return result
    elif Type == 4:
        '''Scipy quadpack'''

        integrand = Integrand(t, r, N, Rho_Liste, w_Liste, K_Liste, kappa, T, Schwartzfunktion)
        tt = time.time()
        Action = integrate.nquad(lambda r1, t1 : integrand(t1,r1),
                [[0,np.pi],[0,T]], opts=[{'epsabs' :10e-10, 'epsrel': 10e-10 },{
                    'epsabs': 10e-10, 'epsrel' : 10e-10} ])

        aa = time.time()
        print('Time it took to integrate in sec:',aa - tt)
        print(Action)
        return Action[0]
    print('done')

def get_integrand_values(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion = True, Comp_String = False, With_Ctypes = True):
    '''
    This method is here for plotting the integrand.
    '''
    if With_Ctypes:
        get_Integrand_with_ctypes(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
            kappa, Schwartzfunktion, Comp_String)
    else:
        get_Integrand_with_c(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String)

    os.system('gcc -o yolo testlib2.c -lm')
    os.system('./yolo')


def Two_Dim_Pic():
    '''
    For N=2 i plotted the Action with respect to K_1 = K_Liste[1] and K_2
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

        norma = get_Action(t, r, N, Integration_bound, T, K_Liste, Rho_Liste,
            w_Liste,kappa, False,False, 1)
        d = open('NumbFor3d.txt', 'w')

        for k1 in K1_Liste:
            for k2 in K2_Liste:
                K_Liste=[k1, k2]
                Wert = get_Action(t, r, N, Integration_bound, T, K_Liste, Rho_Liste, w_Liste,kappa, False,False, 1)
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
    r = si.symarray('r', 1)
    t = si.symarray('t', 1)

    var_K, var_Rho, var_w = configfunktion('Vary_Parameters')
    K_Anf, K_End, K_List = configfunktion('Impuls')
    w_Anf, w_End, w_List = configfunktion('Frequenz')
    Constant, kappa, Rho_List = configfunktion('Constraints')
    Anzahl_N, first = configfunktion('System_sizes')



    T = 1 #Liftime of the universe, it's needed for the Schwartzfunction
    N = 2


    kappa_Anzahl = 1


    x_Anf = 0.
    x_End = np.pi

    Kurve_Names=[]

    Integration_bound = [[x_Anf, x_End], [0, 2*np.pi]]
    Wirk = []
    w_Liste = [1,2,3,4]#eval(w_List)
    K_Liste = [9.7903003749135973, 1.0428185324876262, 0,0.1]
    Rho_Liste = np.array([ 0.23596702, (1- 0.23596702)/3])#,0.1/6,0.1/10 ])#(0.1 +0.03)/6 ])
    Wirkun = get_Action(t, r, N, Integration_bound, T, K_Liste, Rho_Liste,
            w_Liste,kappa, False,False, 1)
    print(Wirkun)

if __name__ == "__main__":
    MainProg()
