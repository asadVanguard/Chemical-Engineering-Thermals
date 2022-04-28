import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sym
from sympy import symbols, Eq, solve
# preparing the class (sympy) in order for it to correspond with the rest of the program
sym.init_printing()

def margules_2_param(x1, x2):
    """
    The goal of this is to figure out the A12 and A21 value, which will then be related to GERT and that will be
    plotted out!
    :param x1: the liquid mole fraction of substance 1
    :param x2: the liquid mole fraction of substance 2
    :param R: The R constant if needed
    :param T: The T value (temperature
    :param PiSat: Saturated pressure that will be resolved from the function- Antoinne's Eq
    :return: returns the Pressure values needed overall for everything else
    """
    A12, A21, GeRTe, GeRTm, P1Sat, P2Sat, gamma1, gamma2 = sym.symbols('A12, A21, GeRTe, GeRTm, '
                                                                       'P1Sat, P2Sat, gamma1, gamma2')

    gamma1 = (x2**2) * (A12 + (2 * (A21 - A12)) * x1)
    gamma2 = (x1**2) * (A21 + (2 * (A12 - A21)) * x2)

    GeRTm = x1 * x2 * ((A12 * x1) + (A21 * x2))

    GeRTe = (x1 * gamma1) + (x2 * gamma2)

    test = solve((GeRTm, GeRTe), (A12, A21, gamma1, gamma2))

    return print(test)


if __name__ == "__main__":
    # code goes here

    x1 = 0.2167
    x2 = 1 - x1

    print(margules_2_param(x1, x2))