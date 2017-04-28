from sympy.printing import ccode
from scipy import integrate
from symengine import I
from SimTest import *
import sympy as sy
import symengine as si
import numpy as np
import random
import time
import os
import Rho_data
import scipy.misc
import ctypes


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
Needed parts for the Integrand
'''

def TensorProduct(mat11, mat22):
    mat1 = np.array(mat11, dtype= object).reshape(2,2)
    mat2 = np.array(mat22, dtype = object).reshape(2,2)
#   print('yooooooooollloo')
#   print(mat1, mat2)
#   print(np.bmat([[mat1[0,0]*mat2, mat1[0,1]*mat2] ,[mat1[1,0]*mat2,
#       mat1[1,1]*mat2]]).reshape(4,4))

    return np.bmat([[mat1[0,0]*mat2, mat1[0,1]*mat2] ,[mat1[1,0]*mat2,
        mat1[1,1]*mat2]]).reshape(4,4)


def prefactor(n):
    return int(scipy.misc.factorial(n + 2))/(8 *
            si.pi**(3/2)*sy.gamma('%d/%d'%(3+2*n,2)))

def diracEigenvalues(n):
    return (2*n + 1)/2

def integralKernelPlus(n, r, theta, phi):
    n=n-1
    lala1 = sy.jacobi(n, 1/2, 3/2, r[0])
    lala2 = sy.jacobi(n, 3/2, 1/2, r[0])
    return prefactor(n)*(si.cos(r[0]/2)*lala1*ident -
                        I*si.sin(r[0]/2)*lala2*sigma_r(theta, phi))

def integralKernelMinus(n, r, theta, phi):
    n=n-1
    lala1 = sy.jacobi(n, 1/2, 3/2, r[0])
    lala2 = sy.jacobi(n, 3/2, 1/2, r[0])
    return prefactor(n)*(si.cos(r[0]/2)*lala1*ident +
                        I*si.sin(r[0]/2)*lala2*sigma_r(theta, phi))

def sigma_r(theta, phi):
    aa = si.sin(theta)*si.cos(phi)*sigma1 + si.sin(theta)*si.sin(phi)*sigma2
    bb =  si.cos(theta)*sigma3

    return aa +bb

def preMatrixPlus(n,K_Liste):
    a = K_Liste[n-1]
    b = si.sqrt(1 + a**2)
    matrix = si.zeros(2)

    matrix[0,0]= 1-b
    matrix[0,1]= a
    matrix[1,0]= -a
    matrix[1,1]= 1+b
    return matrix

def preMatrixMinus(n,K_Liste):
    a = K_Liste[n-1]
    b = si.sqrt(1 + a**2)
    matrix = si.zeros(2)

    matrix[0,0]= 1-b
    matrix[0,1]= -a
    matrix[1,0]= a
    matrix[1,1]= 1+b
    return matrix

def projector(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    mat = np.zeros((4,4), dtype = object)
    for n in range(1,N + 1):
        Koef =Rho_Liste[n-1]*sy.exp(-I*w_Liste[n-1]*t[0])
        #print(Rho_Liste[n-1], sy.exp(-1j*w_Liste[n-1]*t))
        Term1 = TensorProduct(preMatrixPlus(n,K_Liste),integralKernelPlus(n, r, theta, phi))
        Term2 = TensorProduct(preMatrixMinus(n,K_Liste),integralKernelMinus(n, r, theta,phi))
        mat += Koef*(Term1 + Term2)
        #print(Koef, mat)
    return mat

def projectorAdj(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    mat = np.zeros((4,4),  dtype = object)

    for n in range(1, N+1):
        Koeff = Rho_Liste[n-1]*sy.exp(I*w_Liste[n-1]*t[0])
        Term1 = TensorProduct(preMatrixPlus(n,K_Liste),integralKernelMinus(n, r, theta, phi))
        Term2 = TensorProduct(preMatrixMinus(n,K_Liste),integralKernelPlus(n, r,theta, phi))
        mat += Koeff*(Term1 +Term2)
    return mat

def closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste):
    #print(projector(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)*projectorAdj(t, r, theta,phi, N, Rho_Liste, w_Liste, K_Liste))
    return np.dot(projector(t, r, theta, phi, N, Rho_Liste, w_Liste,
        K_Liste),projectorAdj(t, r, theta,phi, N, Rho_Liste, w_Liste, K_Liste))

def lagrangian_without_bound_constr(t, r, theta, phi,N, Rho_Liste, w_Liste, K_Liste):
    sub = closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)
    #print(np.shape(sub),type(sub),sub)
    return np.trace(np.dot(sub,sub)) - 0.25 * np.trace(sub)*np.trace(sub)

'''
Needed parts for the Integrand

'''
'''
Constraints
'''

def boundedness_constraint(t,r,theta, phi, N, Rho_Liste, w_Liste, K_Liste, kappa):
    sub = closedChain(t, r, theta, phi, N, Rho_Liste, w_Liste, K_Liste)
    print ('kappa=', kappa)
    return kappa* np.trace(sub)**2

'''
Constraints
'''

'''
Integrand and Action
'''
def Integrand(t, r, N, Rho_Liste, w_Liste, K_Liste, kappa, T, Schwartzfunktion
        = True):
    args = np.concatenate((t,r))
    exprs = [lagrangian_without_bound_constr(t,r,0,0,N,Rho_Liste,w_Liste,K_Liste)]
    lagr = si.Lambdify(args, exprs, real = False)
    exprs2 = [boundedness_constraint(t,r,0,0,N,Rho_Liste,w_Liste, K_Liste,kappa)]
    bound = si.Lambdify(args,exprs2, real = False)

    if Schwartzfunktion:
        integrand = lambda t1, r1: (max(lagr(t1, r1).real,0) +
                bound(t1,r1).real)*np.sin(r1)**2*np.exp(-(t1)**2/T)
    else:
        integrand = lambda t1, r1 : (max(lagr([t1, r1]).real,0) +
                bound([t1,r1]).real)*np.sin(r1)**2

    return integrand

def get_Integrand_with_c(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String):
    a = ccode(lagrangian_without_bound_constr(t,r,0,0,N, Rho_Liste,
                w_Liste, K_Liste))
    b = ccode(boundedness_constraint(t,r,0,0,N,Rho_Liste,w_Liste, K_Liste,kappa))
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

def get_Test_Integrandt(T):
    Gesamtstring = ''

    Bibliotheken =  '#include <math.h>\n'+'#include <complex.h>\n'+'#include <stdio.h>\n'
    Integranddef = "double f(int n, double args[n])"+ "{return"
    Gesamtstring+= Bibliotheken + Integranddef

    g = open('funcprint.txt', 'r')
    g1 = g.read()

    Func2_1 ='(-1/(cpow(M_PI,0.5)*%2.3f))*'%(T)
    Func2_2 = "cexp(-(cpow(args[0]-1,2)/cpow(%2.3f,2)))"%(T)

    Func2_3 = '+  sin(args[0])'+';'+'}\n'
#    Intdef2 = "double gf(int n, double args[n])"+ "{return 0. ;}"
#    Intdef3 = "double hf(int n, double args[n])"+ "{return 2. ;}"

    Gesamtstring += Func2_1 + Func2_2 + Func2_3 + g1


    f = open('testlib2.c', 'w')
    f.write(Gesamtstring)
    f.close()

    os.system('gcc -x c -shared -o testlib2.so -fPIC testlib2.c')

    g.close()


def get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion = True, Comp_String = False, Type = 1):
    if Type == 1:
        '''Cytpes'''
        get_Integrand_with_ctypes(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
            kappa, Schwartzfunktion, Comp_String)
        aa = time.time()
        lib=ctypes.CDLL('/home/mustafa/Regensburg/Reproduktion_Von_Nikkis_Ergebnissen/Progs_mit_Sympy/testlib2.so')
        lib.f.restype = ctypes.c_double
        lib.f.argtypes = (ctypes.c_int,ctypes.c_double)
        zup = integrate.nquad(lib.f,[[0,np.pi],[0,T]])
        print('(Wirkung, abserr)=',zup)
        handle = lib._handle # obtain the SO handle

        ctypes.cdll.LoadLibrary('libdl.so').dlclose(handle)

        tt = time.time()
        print('Für die Integration benötigte Zeit in sec:',tt-aa)
        return zup
    elif Type ==2:
        '''C'''
        get_Integrand_with_c(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,
        kappa, Schwartzfunktion, Comp_String)
        tt = time.time()
        os.system('./testlib2')
        aa = time.time()
        print(aa-tt)
    elif  Type ==3:
        '''Test'''
        get_Test_Integrandt(T)
        aa = time.time()
        lib=ctypes.CDLL('/home/mustafa/Regensburg/Reproduktion_Von_Nikkis_Ergebnissen/Progs_mit_Sympy/testlib2.so')
        lib.f.restype = ctypes.c_double
        lib.f.argtypes = (ctypes.c_int,ctypes.c_double)
        print(integrate.nquad(lib.f, [[0,np.pi]]))
        tt = time.time()
        print(aa-tt)
    elif Type == 4:
        integrand = Integrand(t, r, N, Rho_Liste, w_Liste, K_Liste, kappa, T, Schwartzfunktion)
        tt = time.time()
        Wirkung = integrate.nquad(lambda r1, t1 : integrand(t1,r1), [[0,np.pi],[0,T]])
        aa = time.time()
        print('Für die Integration benötigte Zeit in sec:',aa - tt)
        print(Wirkung)
        return Wirkung[0]
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

    si.var('theta phi')
    r = si.symarray('r', 1)
    t = si.symarray('t', 1)

    T = 1 #Lebensdauer des Universums, wird fuer die Schwartzfunktion benoetigt
    N = 1

    K_Anzahl=N
    K_Anf = 0.1
    K_End = 1

    kappa = 0.01
    kappa_Anzahl = 1

    w_Liste = [i for i in range(N)]
    Rho_Liste = Rho_data.get_rho_values(N)

    x_Anf = 0.
    x_End = np.pi

    Kurve_Names=[]

    Intgrenze = [x_Anf, x_End]
    Wirk = []
    #get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,kappa, False,False, 4)
    Iter = 5                                      #Die Anzahl der temperaturiterationen
    hilfsarray_fuer_temp = np.linspace(0.01,25,Iter)
    Amplitude = 0.5                                #Amplitude der Fluktuation
                                                   #der Temperatur
    freq = np.pi                                    #Die Frequenz mit der die temp variieren  soll
    halb = 0.001                                      #Der Halbwertswert fuer die
    temp = temperatur(Iter, hilfsarray_fuer_temp, halb, freq, Amplitude)

    '''
    Ich passe das Minimieren zuerst fuer den eindim. Fall an, danach
    passe ich es für mehrere dim. an.
    '''
    x_fitn = random.random()*(K_End - K_Anf) #Startwert... random
    print(x_fitn)
    K_Liste = [x_fitn]
    fitn_wert_x = get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,kappa, False,
            False, 4)
    x_y_min = [x_fitn, fitn_wert_x]
    for m,tt in enumerate(temp):
        for j in range(10):
            randomi = random.random()*(K_End - K_Anf)
            K_Liste = [randomi]
            fitn_wert_y = get_Wirkung(t, r, N, Intgrenze, T, K_Liste, Rho_Liste, w_Liste,     kappa, False, False, 4)
            print('fti_y, fit_y, x = ',fitn_wert_x, fitn_wert_y, x_fitn)
            boltzi = boltzmann(fitn_wert_x, fitn_wert_y, tt)
            print('boltzi = ', boltzi)
            if fitn_wert_x > fitn_wert_y:
                x_fitn = randomi
                x_y_min[0]=x_fitn
                x_y_min[1]=fitn_wert_y
                fitn_wert_x=fitn_wert_y

            elif (fitn_wert_x < fitn_wert_y) and (random.random() < boltzi) :
                x_fitn = randomi
                x_y_min[0]=x_fitn
                x_y_min[1]=fitn_wert_y
                fitn_wert_x = fitn_wert_y
    print(x_y_min[0], x_y_min[1])


