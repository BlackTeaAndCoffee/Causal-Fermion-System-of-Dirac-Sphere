import configparser
import numpy as np
from . import get_config


def configfunktion(section):
    config = configparser.RawConfigParser()
    config.read(get_config('Settings.cfg'))


    if section == 'Vary_Parameters':
        var_rho = config.getboolean('Vary_Parameters', 'Gewichte')
        var_K = config.getboolean('Vary_Parameters', 'Impulse')
        var_w = config.getboolean('Vary_Parameters', 'Frequenz')

        return var_K, var_rho, var_w

    if section == 'Set_Initialstate_randomly':
        rand_Rho = config.getboolean('Set_Initialstate_randomly', 'Weights')
        rand_K =   config.getboolean('Set_Initialstate_randomly', 'Impuls')
        rand_w =   config.getboolean('Set_Initialstate_randomly', 'Frequenz')
        return rand_K, rand_Rho, rand_w

    if section == 'Impuls':

        K_Anf = config.getfloat('Impuls', 'Anf')
        K_End = config.getfloat('Impuls', 'End')
        K_List= config.get('Impuls', 'List')
        return K_Anf, K_End, K_List

    if section == 'Frequenz':

        w_Anf = config.getfloat('Frequenz', 'Anf')
        w_End = config.getfloat('Frequenz', 'End')
        w_List = config.get('Frequenz', 'List')
        return w_Anf, w_End, w_List

    if section == 'Constraints':
        Rho_List = config.get('Constraints', 'List')
        Constant = config.getfloat('Constraints', 'Spur_Konstante')
        kappa = config.getfloat('Constraints', 'Boundary_Konstante')
        return Constant, kappa, Rho_List

    if section == 'System_sizes':

        first = config.getint('System_sizes', 'Erste_Schale')
        Anzahl_N  = config.getint('System_sizes', 'Schalen_Anzahl')

        return Anzahl_N,first
    if section == 'Test':
        Fitness_Nr = config.getint('Test', 'Fitness_Nr')
        StartWithGivenMinima = config.getboolean('Test', 'StartWithGivenMinima')
        return StartWithGivenMinima, Fitness_Nr

