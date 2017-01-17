from sympy import *
from sympy.utilities.lambdify import lambdify
from sympy.abc import r
from sympy.physics.quantum import TensorProduct
from pyx import *
from pyx.graph import axis
from ErstesProgramm import *
import scipy.integrate as integrate
import scipy.special,scipy.misc
import numpy as np
import random
import Rho_data

if __name__ == "__main__":
    N = 2
    Anzahl_an_Ber = 10
    t, theta, phi = symbols('t theta phi')

    K_Anzahl=10
    K_Anf = 0
    K_End = 2.5

    kappa = 0.01
    kappa_Anzahl = 1

    Rho_Liste = Rho_data.get_rho_values(N)
    print(Rho_Liste)

    x_Anf = 0
    x_End = np.pi

    Intgrenze = [x_Anf, x_End]
    Wirk = []

    for i in range(Anzahl_an_Ber):
        w_Liste = np.random.random_sample(N)
        w_Liste[0] = 0
        print(w_Liste)
        K_Liste = np.random.random_sample(N)
        print(K_Liste)
        Wirk.append(get_Wirkung_fuer_kappa(t, r, N, Intgrenze, K_Liste, Rho_Liste, w_Liste,
            kappa))
    print(Wirk)
#   K_Liste =  list(np.linspace(K_Anf,K_End,K_Anzahl))
#   for i in range(K_Anzahl):
#       K2_Liste = [0,K_Liste[i]]
#       Wirk.append(get_Wirkung_fuer_kappa(t, r, N, Intgrenze, K2_Liste, Rho_Liste, w_Liste,
#           kappa))

#   Kurve_Names  = [kappa]
#   diag_plot2(K_Liste, Wirk, 'K', 'Wirkung',
#           Kurve_Names,'N=1_Wirkung_kappa=%0.2f'%kappa, "tr")


