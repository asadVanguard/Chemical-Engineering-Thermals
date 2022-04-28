"""
File:    Equations of States Calculator
Author:  Asad Siddiqui
Date:    10/11/2021
Section: Online Class, 11 AM
E-mail:  asiddiq3@umbc.edu
Description: The goal of this program is to accurately describe the equations of state,
                and condense them into an efficient program
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sym
# preparing the class (sympy) in order for it to correspond with the rest of the program
sym.init_printing()


def van_der_waals(R, temp, press, b, a):
    """
    This function will be modified if not replaced by a more general EOS variable finder, and plotter. For now,
    this is bare bones but will get the job done by far. Eventually, it will take in any EOS, and the function will be
    plotted in goal to find the correct 'independent' variable needed to make it equal to the left side.
    :param R: A universal constant that can be used for any problem
    :param temp: Temperature, in Kelvins
    :param press: Pressure, in the units needed for any problem
    :param b: The specific constant for van der Waals equations (volume of the particles)
    :param a: The specific constant for van der Waals equations (intermolecular attractions on pressure)
    :return: The volume needed to solve the problem (vol) and compressibility factor Z (compressibility_factor)
    """

    # a cubic equation of state for van der Waals defined in terms of Volume
    vol = ((R ** 3 * temp ** 3) / (press * R ** 2 * temp ** 2 + a * press ** 2)) + b

    compressibility_factor = (press * vol) / (R * temp)

    # returning volume calculated
    return vol, compressibility_factor


def truncated_virial_equation(compressibility_factor_z, pressure_p, volume_v, temperature_t, constant_r, func_b, func_c,
                              new_volume_vi, process_code):
    """
    Repetition means extracting the most accurate volume possible, therefore we will extract as needed.
    :param compressibility_factor_z: The compressibility factor Z, if needed be
    :param pressure_p: The pressure of the system, P- only one is really needed to figure out a certain volume
    :param volume_v: The volume of the system V
    :param temperature_t: Temperature T in absolute temperature
    :param constant_r: R constant usually will be 8.314 J/mol*k*m**3
    :param func_b: A function of temperature- can be calculated using this program if needed be
    :param func_c: A function of temperature- can be calculated using this program if needed be
    :param new_volume_vi: The newer volume that will come out and will be reused for the equation
    :param process_code: The process code for which truncated virial equation we need- the Ptizer modified version or no
    :return: The volume desired for this system
    """
    # regular truncation with iteration being done
    if process_code == 0:
        # iteration counter
        i = 0
        # it will iterate up to 3 times
        while i != 3:
            # the equation for regular truncation technique
            new_volume_vi = (constant_r * temperature_t / pressure_p) * (1 + (func_b / volume_v) + (func_c / (volume_v ** 2)))
            # assigning the new term to the older variable to resolve again
            volume_v = new_volume_vi
            # increase the iteration counter by +1
            i += 1
        # your final result!
        new_volume_vi = volume_v

    # Ptizer truncation with iteration being done
    else:
        # we need the Pc to help find Pr and like wise...
        critical_pressure = float(input('What is the critical pressure? '))
        reduced_pressure = pressure_p / critical_pressure
        # Tr will be needed from Tc
        critical_temperature = float(input('What is the critical temperature? '))
        reduced_temperature = temperature_t / critical_temperature
        # omega will be needed to help calculate B-hat
        omega = float(input('What is the value of the Accentric Value? '))

        # calculating B-naught
        b_naught = 0.083 - (0.422 / (reduced_temperature ** 1.6))
        # calculating B-one
        b_one = 0.139 - (0.172 / (reduced_temperature ** 4.2))
        # calculating factor Z using the equation from the truncated equation
        compressibility_factor_z = 1 + (b_naught * (reduced_pressure / reduced_temperature)) + (b_one * (reduced_pressure/ reduced_temperature))

        # calculating the newer volume that the user needs at last by rearranging the ideal gas Z to truncated Z
        new_volume_vi = ((constant_r * temperature_t) * compressibility_factor_z) / temperature_t

    # return the final answer for the player
    return compressibility_factor_z, new_volume_vi


def redlich_kwong(compressibility_factor_z, pressure_p, volume_v, temperature_t, constant_r, func_b, func_a):
    """
    The Redlich-Kwong equation.
    :param compressibility_factor_z: The compressibility factor Z, if needed be
    :param pressure_p: The pressure of the system, P- only one is really needed to figure out a certain volume
    :param volume_v: The volume of the system V
    :param temperature_t: Temperature T in absolute temperature
    :param constant_r: R constant usually will be 8.314 J/mol*k*m**3
    :param func_b: A function of temperature- can be calculated using this program if needed be
    :param func_a: A function of temperature- can be calculated using this program if needed be
    :return: The volume and Compressibility Factor Z values
    """
    # asking user for the critical temperature input of their choice
    crit_temp = float(input('Furthermore, what is the value of the critical temperature? '))
    # asking user for the critical pressure input of their choice
    crit_press = float(input('Lastly, what is the value of the critical pressure? '))
    # another methodology for calculating and solving equations, via SymPy. Very advanced tool.
    Z, P, V, T, R, Tc, Pc, func_a, func_b, A, B = sym.symbols('Z, P, V, T, R, Tc, Pc, func_a, func_b, A, B')
    # assigning several variables their respective values
    P = pressure_p
    T = temperature_t
    R = constant_r
    Tc = crit_temp
    Pc = crit_press
    # solving for a and b, since they are products of constants, and Pc and Tc
    func_a = 0.42798 * (((R ** 2) * (Tc ** 2.5)) / Pc)
    func_b = 0.08664 * (R * Tc / Pc)
    # solving for A and B since they are functions of temperature and help get us out respective values for Z
    A = (func_a * P) / ((R ** 2) * (T ** 1.5))
    B = (func_b * P) / ((R) * (T))

    # the methodology for SymPy asks us to essentially create the equation, and declare what must be solved... volume!
    volume = sym.Eq(P - (((R * T) / (V - func_b)) -
                                   (func_a / ((T ** 0.5) * V * (V + func_b)))))
    # solving for Z factor with the same method
    compressibility_factor_z = sym.Eq((Z ** 3) - (Z ** 2) + (A - B - (B ** 2)) * Z - (A * B))

    # the method inherently is recursive, and so we only want the first 1 or 2 recursions of said functions to return
    return print(sym.solve(volume, V)), print(sym.solve(compressibility_factor_z, Z))


def soave_redlich_kwong(pressure_p, temperature_t, constant_r):
    """
    The Soave-Redlich-Kwong equation.
    :param pressure_p: The pressure of the system, P- only one is really needed to figure out a certain volume
    :param temperature_t: Temperature T in absolute temperature
    :param constant_r: R constant usually will be 8.314 J/mol*k*m**3
    :return: The volume and Compressibility Factor Z values
    """
    # another methodology for calculating and solving equations, via SymPy. Alpha, Omega, and m are needed to resolve.
    Z, P, V, T, R, Tc, Pc, func_a, func_b, A, B = sym.symbols('Z, P, V, T, R, Tc, Pc, func_a, func_b')
    # asking user for the critical temperature input of their choice
    crit_temp = float(input('Furthermore, what is the value of the critical temperature? '))
    # asking user for the critical pressure input of their choice
    crit_press = float(input('Next, what is the value of the critical pressure? '))
    # asking user for an input of accentric value, omega
    omega = float(input('Lastly, what is the accentric value? (OMEGA) '))

    # solving for reduced temperature to help out with future calculations
    Tr = temperature_t / crit_temp
    # solving for reduced pressure to help out with future calculations
    Pr = pressure_p / crit_press
    # solving for variable-m that will then be used to solve for alpha! Takes into consideration of accentric values
    m = 0.48508 + (1.579 * omega) - ((0.176 * omega) ** 2)
    # solving for alpha which will be used for calculations for function-a and other calculations
    alpha = (1 + m * (1 - np.sqrt(Tr))) ** 2

    # assigning several variables their respective values
    P = pressure_p
    T = temperature_t
    R = constant_r
    Tc = crit_temp
    Pc = crit_press
    W = omega
    Al = alpha
    M = m

    # solving for a and b, since they are products of constants, and Pc and Tc
    func_a = ((0.42798 * Al * ((R * Tc) ** 2)) / Pc)
    func_b = (0.08664 * R * Tc) / Pc

    # solving for A and B since they are functions of temperature and help get us out respective values for Z
    A = (func_a * P) / ((R * T) ** 2)
    B = (func_b * P) / (R * T)

    # the methodology for SymPy asks us to essentially create the equation, and declare what must be solved... volume!
    volume = sym.Eq(P - (((R * T) / (V - func_b)) -
                         (func_a / (V * (V + func_b)))))
    # solving for Z factor with the same method
    compressibility_factor_z = sym.Eq((Z ** 3) - (Z ** 2) + (A - B - (B ** 2)) * Z - (A * B))

    # the method inherently is recursive, and so we only want the first 1 or 2 recursions of said functions to return
    return print(sym.solve(volume, V)), print(sym.solve(compressibility_factor_z, Z))


def peng_robinson(pressure_p, temperature_t, constant_r):
    """
    The Peng Robinson Equation of State
    :param pressure_p: The pressure of the system
    :param temperature_t: The temperature of the system
    :param constant_r: The R-constant
    :return: Returns the needed Z-factor and Volume
    """
    # defining the symbols for sympy to recognize as we solve the equations
    V, Z, P, R, T, a, b = sym.symbols('V, Z, P, R, T, a, b')
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

    # assigning several variables that will finally be used for the equation solves to their distinctive values
    P = pressure_p
    T = temperature_t
    R = constant_r
    a = func_a
    b = func_b

    # the equation solver for finding volume can be done with sym.Eq-- the equation is Peng/Robinson
    volume = sym.Eq(P - (((R * T) / (V - b)) - (a / ((V * (V + b)) + (b * (V - b))))))

    volume2 = sym.solve(volume, V)

    # the equation solver for finding the compressibility factor Z, except this will return in terms of V
    compressibility_z = sym.Eq(Z - (volume2[0]/ (volume2[0] - b)) - ((a * volume2[0]) /(R * T)) * (1 / ((volume2[0] * (volume2[0] + b)) + (volume2[0] * (volume2[0] - b)))))

    # return the user the resolved values of volume compressibility factor Z; only the first recursion is needed!
    return print(sym.solve(volume, V)), print(sym.solve(compressibility_z, Z))


# main function
if __name__ == '__main__':
    # asking user for input if they would like to use the program
    command_center = input('Would you like to calculate the volume from your state function of choice? (yes/no) ')

    # run conditional one if 'yes'
    if command_center == 'yes':
        command_center = input('What is the equation you want to solve for? Scroll up for a list of options: ')

        if command_center == 'van der Waals':
            # variable to be used for cutting off while loop to prevent infinite looping
            restart_prompt = input('What is the value of R, or type cancel to back out right now: ')
            # creating an empty list that will be filled up if needed be, or will remain harmless or empty
            list_of_terms = []

            # run while loop until user says 'cancel'
            while restart_prompt != 'cancel':
                # defining whatever the user input into a decimal point number
                R = float(restart_prompt)
                # the following are simply strings for each variable to help distinguish which one is which later on
                temp = 'Temperature'
                press = 'Pressure'
                b = 'b-constant'
                a = 'a-constant'

                # creating a list for the variables that will be used for a van der Waals volume finder
                list_of_terms = [R, temp, press, b, a]

                # running a for-loop just from the start of 'temp' through 'a'. R is not needed since it was
                # already asked..
                for i in range(1, len(list_of_terms)):
                    # asking user for number input for each variable needed in accordance from the list- convert
                    # to decimal
                    def_variables = float(input('What is the value of {}? '.format(list_of_terms[i])))

                    # readjusting the strings that belonged to each variable to the decimal number needed-
                    # tiny inefficient
                    list_of_terms[i] = def_variables

                # asking user if they would like to continue, or stop. Either way, it will tell the volume immediately
                restart_prompt = input('Would you like to continue, or cancel? Enter cancel to stop, '
                                       'or input another R value. ')

                # giving the user their volume desired
                print(van_der_waals(list_of_terms[0], list_of_terms[1], list_of_terms[2], list_of_terms[3],
                                    list_of_terms[4]))

        elif command_center == 'truncated virial equation':
            order = input('Do you have B and C already? That is another format of the EOS. ')

            if order == 'yes':
                # variable to be used for cutting off loop to prevent infinite looping
                new_prompt = float(input('What is the value of R, or type cancel to back out right now: '))
                # list of variables to set up out equation of state
                compressibility_z = 'Z'
                press_p = 'P'
                vol_v = 'V'
                temp_t = 'T'
                const_r = 'R'
                func_b = 'B'
                func_c = 'C'
                new_vol_vi = 'Vi'
                # process code to distinguish what we need
                process_code = 'Enter the number for the process code- 0 for regular truncation.'

                # creating an empty list that will be filled up if needed be, or will remain harmless or empty
                list_of_terms = [compressibility_z, press_p, vol_v, temp_t, const_r, func_b, func_c, new_vol_vi,
                                 process_code]

                # creating a for-loop to set up the exact numbers needed for this equation of state
                for i in range(len(list_of_terms)):
                    # assigning decimal numbers to each variable
                    def_variables = float(input('Alright, input a value for... {}! '.format(list_of_terms[i])))
                    # reassigning our list with newer terms
                    list_of_terms[i] = def_variables
                # evaluating the new volume for the user to finally use
                print(truncated_virial_equation(list_of_terms[0], list_of_terms[1], list_of_terms[2], list_of_terms[3],
                                                list_of_terms[4], list_of_terms[5], list_of_terms[6], list_of_terms[7],
                                                list_of_terms[8]))

            # otherwise, this is process code 1 whereas we need accentric values, and other stuff
            else:
                # a list of variables to set up our EOS
                compressibility_z = 'Z - enter 0'
                press_p = 'P'
                vol_v = 'V'
                temp_t = 'T'
                const_r = 'R'
                func_b = 'B, enter 0'
                func_c = 'C, enter 0'
                new_vol_vi = 'Vi - 0'
                # process code for our truncated virial equation with accentric values
                process_code = 'Enter the number for the process code: 1 for Ptizer mod. '

                # creating an empty list that will be filled up if needed be, or will remain harmless or empty
                list_of_terms = [compressibility_z, press_p, vol_v, temp_t, const_r, func_b, func_c, new_vol_vi,
                                 process_code]

                # creating a for-loop to set up the exact numbers needed for this equation of state
                for i in range(len(list_of_terms)):
                    # assigning decimal numbers to each variable
                    def_variables = float(input('Alright, input a value for... {}! '.format(list_of_terms[i])))
                    # reassigning our list with newer terms
                    list_of_terms[i] = def_variables

                # evaluating the new volume for the user to finally use
                print(truncated_virial_equation(list_of_terms[0], list_of_terms[1], list_of_terms[2], list_of_terms[3],
                                                list_of_terms[4], list_of_terms[5], list_of_terms[6], list_of_terms[7],
                                                list_of_terms[8]))
        elif command_center == 'redlich kwong':
            # list of variables to set up out equation of state
            compressibility_z = 'Z -  enter 0'
            pressure_p = 'P'
            vol_v_guess = 'V - enter 0'
            temperature_t = 'T'
            constant_r = 'R'
            func_b = 'b- enter 0'
            func_a = 'a- enter 0'

            # creating an empty list that will be filled up if needed be, or will remain harmless or empty
            list_of_terms = [compressibility_z, pressure_p, vol_v_guess, temperature_t, constant_r, func_b, func_a]
            # creating a for-loop to set up the exact numbers needed for this equation of state
            for i in range(len(list_of_terms)):
                # assigning decimal numbers to each variable
                def_variables = float(input('Alright, input a value for... {}. '.format(list_of_terms[i])))
                # reassigning our list with newer terms
                list_of_terms[i] = def_variables
            # print out the result of the Redlich-Kwong equation of state-- the volume and compressibility factor!
            print(redlich_kwong(list_of_terms[0], list_of_terms[1], list_of_terms[2], list_of_terms[3], list_of_terms[4],
                                list_of_terms[5], list_of_terms[6]))

        elif command_center == 'SRK':
            # list of variables to set up out equation of state
            pressure_p = 'P'
            temperature_t = 'T'
            constant_r = 'R'

            # creating an empty list that will be filled up
            list_of_terms = [pressure_p, temperature_t, constant_r]
            # creating a for-loop to set up the exact numbers needed for this equation of state
            for i in range(len(list_of_terms)):
                # assigning decimal numbers to each variable
                def_variables = float(input('Alright, input a value for... {}. '.format(list_of_terms[i])))
                # reassigning our list with newer terms
                list_of_terms[i] = def_variables
            # print out the result of the Soave-Redlich-Kwong equation of state-- volume and compressibility factor!
            print(soave_redlich_kwong(list_of_terms[0], list_of_terms[1], list_of_terms[2]))

        elif command_center == 'P/R':
            # list of variables to set up our equation of state
            pressure_p = 'P'
            temperature_t = 'T'
            constant_r = 'R'

            # creating an empty list that will be filled up
            list_of_terms = [pressure_p, temperature_t, constant_r]
            # creating a for-loop to set up the exact numbers needed for this EOS
            for i in range(len(list_of_terms)):
                # assigning decimal numbers to each variable
                def_variables = float(input('Alright, input a value for... {}. '.format(list_of_terms[i])))
                # reassigning our list with newer terms
                list_of_terms[i] = def_variables
            # printing out the results of the Peng and Robinson Equation of State-- i.e, the volume and compressibility!
            print(peng_robinson(list_of_terms[0], list_of_terms[1], list_of_terms[2]))