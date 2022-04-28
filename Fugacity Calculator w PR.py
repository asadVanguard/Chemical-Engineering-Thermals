"""
File:    Fugacity Calculator w/ Peng & Robinson
Author:  Asad Siddiqui
Date:    10/28/2021
Section: Online Class, 11 AM
E-mail:  asiddiq3@umbc.edu
Description: The goal of this calculator is to essentially calculate the fugacity value using the EOS of Peng & Robinson
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sym
# preparing the class (sympy) in order for it to correspond with the rest of the program
sym.init_printing()


def peng_robinson(pressure_p, temperature_t, constant_r):
    """
    The Peng Robinson Equation of State
    :param pressure_p: The pressure of the system
    :param temperature_t: The temperature of the system
    :param constant_r: The R-constant
    :return: Returns the needed Z-factor and Volume
    """
    # defining the symbols for sympy to recognize as we solve the equations
    V, Z, P, R, T, a, b, A, B, f = sym.symbols('V, Z, P, R, T, a, b, A, B, f')
    # asking user for the critical temperature
    crit_temp = float(input('Alright, what is the critical temperature needed? '))
    # asking user for the ciritcal pressure
    crit_press = float(input('Next, what is the critical pressure needed? '))
    # asking user for the accentric value
    omega = float(input('Lastly, what is the accentric value for omega? '))
    # finding the reduced temperature for future calculations
    Tr = temperature_t / crit_temp
    # defining variable 'b' that is related to the equation
    func_b = (0.07780 * constant_r * crit_temp) / crit_press
    # defining variable 'kappa' and using it for future calculations via the accentric value
    kappa = 0.37469 + (1.54226 * omega) - (0.2669 * (omega ** 2))
    # defining alpha via the use of kappa and reduced temperature
    alpha = (1 + (kappa * (1 - (Tr ** 0.5))))**2
    # solving for the critical value of our function of temperature, 'a'
    crit_a = 0.45724 * ((constant_r ** 2) * (crit_temp ** 2)/ crit_press)
    # finally, solving for function of temperature 'a'
    func_a = alpha * crit_a

    cap_a = (func_a * pressure_p) / (constant_r**2 * temperature_t**2 )

    cap_b = (func_b * pressure_p) / (constant_r * temperature_t)
    # assigning several variables that will finally be used for the equation solves to their distinctive values
    P = pressure_p
    T = temperature_t
    R = constant_r
    a = func_a
    b = func_b
    A = cap_a
    B = cap_b
    # the equation solver for finding volume can be done with sym.Eq-- the equation is Peng/Robinson
    volume = sym.Eq(P - (((R * T) / (V - b)) - (a / ((V * (V + b)) + (b * (V - b))))))

    volume2 = sym.solve(volume, V)

    # the equation solver for finding the compressibility factor Z, except this will return in terms of V
    compressibility_z = sym.Eq(Z - (volume2[0]/ (volume2[0] - b)) -
                               ((a * volume2[0]) /(R * T)) * (1 / ((volume2[0] * (volume2[0] + b)) +
                                                                   (volume2[0] * (volume2[0] - b)))))

    # solving for the compressibility factor and storing it into a dummy variable
    compressibility_z_2 = sym.solve(compressibility_z, Z)

    # removing the tailing imaginary number since it is negligible, and subtracting it by 1 to ideally get around
    # the same compressibility factor that is desired.
    compressibility_z_3 = 1 - abs((compressibility_z_2[0] * (10 ** -2)))

    # this is the fugacity equation with respect to Peng & Robinson
    fugacity = P * math.exp((compressibility_z_3 - 1) -
                            (sym.log(compressibility_z_3 - cap_b) - (cap_a / (2 * cap_b * (2 ** (1 / 2))))) *
                            sym.log((compressibility_z_3 + (1 + 2 ** (1 / 2)) * cap_b) / (compressibility_z_3 +
                                                                                          (1 - 2 ** (1 / 2)) * cap_b)))
    return fugacity


# main function
if __name__ == '__main__':

    # list of variables to set up our equation of state and our list (inefficient, but gets the job done)
    pressure_p = 'P'
    temperature_t = 'T'
    constant_r = 'R'
    point_1 = ''
    point_2 = ''
    point_3 = ''
    point_4 = ''
    point_5 = ''
    point_6 = ''
    point_7 = ''
    point_8 = ''
    point_9 = ''
    point_10 = ''

    # creating an empty list that will be filled up
    list_of_terms = [pressure_p, temperature_t, constant_r]
    list_of_fugacities = [point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8, point_9, point_10]
    list_of_pressures = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # creating an iterator for the while loop to run off of
    i = 0
    while i != 10:
        # the terms needed for Peng/Robinson are going to be distributed among multiple pressure values
        # that we are using
        list_of_terms[0] = list_of_pressures[i]
        # this is the temperature of the problem
        list_of_terms[1] = 373.15
        # this is the R constant of the problem
        list_of_terms[2] = 83.14
        # creating a list of fugacity terms depending on what the peng_robinson equation comes out with
        list_of_fugacities[i] = peng_robinson(list_of_terms[0], list_of_terms[1], list_of_terms[2])
        # increasing our iterator counter by +1 each time
        i += 1

    # a little print out of what comes out in the end and is being placed into the graph itself
    print(list_of_fugacities, list_of_pressures)

    """
    The following below are all a part of PLT, or MATPLOT import into Python. Essentially we are plotting our lists!
    """
    plt.plot(list_of_pressures, list_of_fugacities, r'^')
    plt.xlabel("Pressure (Bars)")
    plt.ylabel("Fugacities (Bars)")
    plt.title("Peng & Robinson Fugacity v. Pressure")
    plt.xlim(5, 110, 10)
    plt.ylim(5, 265, 10)
    plt.show()