import random
import numpy as np


def boltzmann(f_x, f_y, temperatur, Boltzmann_Konstant):
    print('f_x,f_y', f_x, f_y)
    diff = abs(f_y - f_x)
    print('diff', diff)
    wert = np.exp(- diff/(temperatur*Boltzmann_Konstant))
    return wert

def Temperatur_Function(tempus, decay_const, freq, Amplitude):
    return np.exp(-decay_const*tempus**2)*(Amplitude*np.cos(freq*tempus*0.5) + 1)

def temperatur(Iter, hilfsarray_fuer_temp, decay_const, freq, Amplitude):
    temperatur_list =[]
    for k in range(Iter):
        temperatur_list.append(Temperatur_Function(hilfsarray_fuer_temp[k],decay_const,freq,
            Amplitude))
    return temperatur_list

