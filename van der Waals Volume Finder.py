"""
File:    van der Waals Volume Finder
Author:  Asad Siddiqui
Date:    10/11/2021
Section: Online Class, 11 AM
E-mail:  asiddiq3@umbc.edu
Description: The goal of this program is to solve for the volume from this cubic state of function via Python
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


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


def redlich_kwong(volume_v):
    """
    The Redlich-Kwong equation.
    :param compressibility_factor_z: The compressibility factor Z, if needed be
    :param pressure_p: The pressure of the system, P- only one is really needed to figure out a certain volume
    :param volume_v: The volume of the system V (guess to make the solver work)
    :param temperature_t: Temperature T in absolute temperature
    :param constant_r: R constant usually will be 8.314 J/mol*k*m**3
    :param func_b: A function of temperature- can be calculated using this program if needed be
    :param func_a: A function of temperature- can be calculated using this program if needed be
    :return: The volume and Compressibility Factor Z values
    """

    return list_of_terms[1] - ((list_of_terms[4] * list_of_terms[3]) / (volume_v - list_of_terms[5])) - \
           (list_of_terms[6] / ((list_of_terms[3] ** 0.5) * volume_v * (volume_v + list_of_terms[5])))




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
                # variable to be used for cutting off while loop to prevent infinite looping
                new_prompt = float(input('What is the value of R, or type cancel to back out right now: '))
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
            compressibility_z = 'Z'
            pressure_p = 'P'
            vol_v_guess = 'V - enter 1'
            temperature_t = 'T'
            constant_r = 'R'
            func_b = 'b'
            func_a = 'a'

            # creating an empty list that will be filled up if needed be, or will remain harmless or empty
            list_of_terms = [compressibility_z, pressure_p, vol_v_guess, temperature_t, constant_r, func_b, func_a]
            # creating a for-loop to set up the exact numbers needed for this equation of state
            for i in range(len(list_of_terms)):
                # assigning decimal numbers to each variable
                def_variables = float(input('Alright, input a value for... {}. '.format(list_of_terms[i])))
                # reassigning our list with newer terms
                list_of_terms[i] = def_variables

            critical_temperature = float(input('What is the value of the critical temperature? '))

            critical_pressure = float(input('What is the value of the critical pressure? '))

            func_a = 0.42798 * (((list_of_terms[4] ** 2) * (critical_temperature ** 2.5)) / critical_pressure)

            func_b = 0.08664 * (list_of_terms[4] * critical_temperature / critical_pressure)

            volume_v = np.linspace(0, 0.005)

            volume_v = fsolve(redlich_kwong, volume_v)

            print(volume_v)