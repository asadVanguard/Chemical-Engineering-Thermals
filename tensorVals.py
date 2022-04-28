"""
File:   tensorVals.py
Author:  Asad Siddiqui
Date:    3/25/20222
Section: Section 2, Tuesday/Thursday 11:30 AM
E-mail:  asiddiq3@umbc.edu
Description: This file is meant for all calculations pertaining to convection heat transfer. Use at discretion.
"""

import openpyxl
from openpyxl import load_workbook


# function for linearly interpolating data that might not be catered to one's needs
def lin_interpolate(x, x1, x2, y1, y2):
    """
    The goal of this program is to take nay set of values needed, from a PROPERLY formatted Excel Spreadsheet
    :param x: The desired temperature that will have its new interpolated value
    :param x1: The lower temperature
    :param x2: The higher temperature
    :param y1: The lower temperature's respective Appendix Piece
    :param y2: The higher temperature's respective Appendix Piece
    :return: The interpolated value (y) for what is needed
    """
    vals_needed = input("What values do you need interpolated in particular, or say \'all' if you need every value"
                        "interpolated? (One at a time, please) If you need nothing else, say \'stop.': ")
    if vals_needed == "all":
        for i in range(len(y1)):
            y_inter = y1[i] + (x - x1) * ((y2[i] - y1[i]) / (x2 - x1))
            interp_value.append(y_inter)
        return interp_value
    elif vals_needed == "density":
        interp_value.append(y1[1] + (x - x1) * ((y2[1] - y1[1]) / (x2 - x1)))
        return lin_interpolate(tfilm, temp_1, temp_2, y1, y2)
    elif vals_needed == "specific heat":
        interp_value.append(y1[2] + (x - x1) * ((y2[2] - y1[2]) / (x2 - x1)))
        return lin_interpolate(tfilm, temp_1, temp_2, y1, y2)
    elif vals_needed == "thermal conductivity":
        interp_value.append(y1[3] + (x - x1) * ((y2[3] - y1[3]) / (x2 - x1)))
        return lin_interpolate(tfilm, temp_1, temp_2, y1, y2)
    elif vals_needed == "thermal diffusivity":
        interp_value.append(y1[4] + (x - x1) * ((y2[4] - y1[4]) / (x2 - x1)))
        return lin_interpolate(tfilm, temp_1, temp_2, y1, y2)
    elif vals_needed == "dynamic_viscosity":
        interp_value.append(y1[5] + (x - x1) * ((y2[5] - y1[5]) / (x2 - x1)))
        return lin_interpolate(tfilm, temp_1, temp_2, y1, y2)
    elif vals_needed == "kinematic viscosity":
        interp_value.append(y1[6] + (x - x1) * ((y2[6] - y1[6]) / (x2 - x1)))
        return lin_interpolate(tfilm, temp_1, temp_2, y1, y2)
    elif vals_needed == "prandtl number":
        interp_value.append(y1[7] + (x - x1) * ((y2[7] - y1[7]) / (x2 - x1)))
        return lin_interpolate(tfilm, temp_1, temp_2, y1, y2)
    else:
        return interp_value


# function for finding the value of the coefficient of convection, Ra, and Nu
def natural_conv_h(L, W, Diameter, theta, k, g, beta, Pr, kin_visc, Ts, Tconvec):
    """
    Function is meant to find the value of "h" when it is natural convection dominant, and will be referred to for the
    calculations of h-natural, and Ra, as well as Nusselts Number
    :param L: Length of an object, preferably a rectangular surface or cylindrical surface
    :param W: Width of an object, preferably a rectangular surface
    :param Diameter: Diameter of an object, preferably a sphere or cylindrical object
    :param theta: The angle at which the object has been tilted, and must be resolved with our natural convection flow
    :param k: coefficient of conduction
    :param g: gravitational constant
    :param beta: coefficient of volume expansion, or 1/K, or essentially, 1/T if we assume ideal gas
    :param Pr: Prandtl's Number
    :param kin_visc: kinematic viscosity of a fluid
    :param Ts: Surface temperature
    :param Tconvec: Temperature of the fluid moving across the surface of the object's Ts (the desired Ts of focus)
    :return: the value of "h" for further calculations.
    """

    # asking user for the geometry of the shape needed
    geometry_question = input("Describe the shape: a vertical plate, inclined plate, horizontal plate, "
                              "vertical cylinder, horizontal cylinder, or a sphere. (LOWER CASE, INCLUDE SPACE): ")
    # creating a series of if statements to help find the best Ra equation needed
    if geometry_question == "vertical plate":
        # setting our characteristic length to be equal to length in this scenario
        Lc = L
        # that said Ra equation is the following below:
        rayleigh_num = abs((g * beta * (Ts - Tconvec) * (Lc ** 3) / (kin_visc ** 2)) * Pr)
        # creating another set of if statements to pinpoint what nusselts number is exactly needed
        if rayleigh_num >= 1*10**4 or rayleigh_num <= 1*10**9:
            nusselts_num = 0.59*(rayleigh_num**(1/4))
        else:
            nusselts_num = 0.1*(rayleigh_num**(1/3))
    # horizontal plate is the other condition, at this current point of time. Will NOT be necessary for the time being
    elif geometry_question == "horizontal plate":
        # the characteristic length of a horizontal plate scenario
        Lc = (L * W) / (2*L + 2*W)
        # the respective Ra equation needed for such a scenario
        rayleigh_num = (g * beta * (Ts - Tconvec) * (Lc ** 3) / (kin_visc ** 2)) * Pr
        # another series of if-statements to help pinpoint what Nusselts Number eq. is needed
        if rayleigh_num >= 1*10**4 or rayleigh_num < 1*10**7:
            nusselts_num = 0.59*(rayleigh_num**(1/4))
        elif rayleigh_num >= 1*10**7 or rayleigh_num <= 1*10**11:
            nusselts_num = 0.1*(rayleigh_num**(1/3))
    elif geometry_question == "horizontal cylinder":
        Lc = Diameter / 2
        rayleigh_num = (g * beta * (Ts - Tconvec) * (Lc ** 3) / (kin_visc ** 2)) * Pr

        nusselts_num = (0.6 + (0.387*rayleigh_num**(1/6)) / ((1 + (0.559/Pr)**(9/16))**(8/27)))**2
    # after finding the value of Nusselts Number, we will now find the coefficient of natural convection
    h = (nusselts_num * k) / L
    # return a tuple of coefficient of convection, and Ra
    return h, rayleigh_num, nusselts_num


# function meant to find the value of the coefficient of convection for window scenarios, and more enclosed scenarios
def enclosure_checker(L, H, W, k, g, Pr, kin_visc, beta, Tsin, Tsout):
    """
    This will optimize our enclosed space, and help us find the coefficient of convection needed for such scenarios-
    this function is ONLY for the rectangular variants at the moment.
    :param L: length of the enclosed figure
    :param H: height of the enclosed figure
    :param W: width of the enclosed figure
    :param k: coefficient of conduction
    :param Pr: Prandtl's Number
    :param kin_visc: kinematic viscosity of a fluid
    :param beta: coefficient of volume expansion, or 1/K, or essentially, 1/T if we assume ideal gas
    :param Tsin: Inner surface temperature
    :param Tsout: Outer surface temperature
    :return: The value of the coefficient of convection
    """
    # asking the user for the type of enclosed surface it is
    type_of_enclosure = input("What is the type of enclosed space that you want to optimize the coefficient of "
                              "convection for? (horizontal, inclined, or vertical?) ")
    # setting our characteristic length to be equal to length, L in this case since we will only examine 3 rec. types.
    Lc = L
    # the respective Ra number equation that will be used for this scenario
    rayleigh_num = abs((g * beta * (Tsin - Tsout) * (Lc ** 3) / (kin_visc**2)) * Pr)
    # setting up a series of if-statements to pinpoint the exact Nusselts Number equation needed
    if type_of_enclosure == "vertical" or type_of_enclosure == "vertical rectangular enclosure":
        # 3 conditions will be used to find what the proper equation is. In this case, the first is:
        cond_1 = H/L
        # the second condition is Prandtl's Number
        cond_2 = Pr
        # the last condition is Ra's number
        cond_3 = rayleigh_num
        # using the parameters specified by Cengal, et. al. We can make another series of nested if-statements
        if 2 > cond_1 >= 1:
            # for this case, the Nu equation shall be this:
            nusselts_num = 0.18 * (((Pr / (0.2 + Pr)) * rayleigh_num)**0.29)
        # case 2 specified by Cengal, et. al.
        elif (10 > cond_1 >= 2) or cond_3 < (1 * 10 ** 10):
            # this will be the equation of Nu for case 2
            nusselts_num = 0.22 * (((Pr / (0.2 + Pr)) * rayleigh_num)**0.28) * ((H/L)**(-1/4))
        # case 3 specified
        elif (40 > cond_1 >= 10) or ((1 * 10 ** 7) > cond_3 >= (1 * 10 ** 4)):
            # this will be the equation of Nu for case 3
            nusselts_num = 0.42 * (rayleigh_num**(1/4)) * (Pr**0.012) * ((H/L)**(-0.3))
        # otherwise, this will be the last case scenario
        else:
            # the final possible equation
            nusselts_num = 0.046 * rayleigh_num**(1/3)
    # if the user made an error, they will be brought up with this message.
    else:
        print("try again")
        # return a tuple full of values that can be reused for the user
        return enclosure_checker(L, H, W, k, g, Pr, kin_visc, beta, Tsin, Tsout)
    # if all else was successful, we will now find the desired coefficient of convection for our enclosed surface
    h = (k * nusselts_num / L)
    # and now, we will return a tuple for the user to use for future reasons
    return h, cond_1, cond_3


# main block of code
if __name__ == "__main__":
    # an empty list that can be appended for future uses
    interp_value = []
    # extracting an Excel Spreadsheet that has access to the data in Cengal, et. al. Appendices
    tensor_book = load_workbook("E:\\Excel and Py\\TensorVals.xlsx")
    # letting Python read the Excel spreadsheet now
    tensor_page = tensor_book.active
    # welcoming our user for the program
    print('Welcome to the Temperature Sensor Value program, where our goal is to simply deliver the interpolated values'
          'of the data for any gas species/air input into our system.')
    # asking user for the range of data that must be accessed from the Excel Spreadsheet
    data_range1 = input("What is the lower limit row of data needed for this? (List the # only please): ")
    data_range2 = input("How about the upper limit row of data? (List the # only please): ")
    # asking the user if they need to input temperatures that are not the default ones used in the textbook
    mod_temps = input("Even though you selected a data range, is one of your temps actually not the same as the data"
                      "range's temperature scales? (y/n) ")
    # string concatenation
    tensor_vals = tensor_page['A' + data_range1:'H' + data_range1]
    tensor_vals2 = tensor_page['A' + data_range2:'H' + data_range2]
    # creating a for-each style loop to extract each value in each data-cell from Excel
    for c1, c2, c3, c4, c5, c6, c7, c8 in tensor_vals:
        temp_1 = c1.value
        dens_1 = c2.value
        specheat_1 = c3.value
        thermalcond_1 = c4.value
        thermaldiffus_1 = c5.value
        dynamicvisc_1 = c6.value
        kinematicvisc_1 = c7.value
        prandtl_1 = c8.value
    # once again, doing the same but for the higher range of values
    for c1, c2, c3, c4, c5, c6, c7, c8 in tensor_vals2:
        temp_2 = c1.value
        dens_2 = c2.value
        specheat_2 = c3.value
        thermalcond_2 = c4.value
        thermaldiffus_2 = c5.value
        dynamicvisc_2 = c6.value
        kinematicvisc_2 = c7.value
        prandtl_2 = c8.value
    # if the user DID need to input newer temperature values for the process,
    if mod_temps == 'y':
        # we will now ask them for the lower temperature in the case
        temp_1 = float(input("What is the lower value temperature? "))
        # and ask for the higher temperature.
        temp_2 = float(input("What is the upper value temperature? "))
    # if they DID NOT
    else:
        # we will now continue with the rest of the process
        print("Continuing with the process now... hold on a moment.")
    # calculating the average temperature in contact with the air for future convection calculations
    tfilm = (temp_1 + temp_2)/2
    # creating 2 lists that has the values extracted from the Excel sheet used
    y1 = [dens_1, specheat_1,thermalcond_1, thermaldiffus_1, dynamicvisc_1, kinematicvisc_1, prandtl_1]
    y2 = [dens_2, specheat_2, thermalcond_2, thermaldiffus_2, dynamicvisc_2, kinematicvisc_2, prandtl_2]

    # CALLING A FUNCTION TO NOW INTERPOLATE THE DATA FOR MUCH MORE ACCURACY
    values_returned = lin_interpolate(tfilm, temp_1, temp_2, y1, y2)

    # asking the user for what they would like to do next
    further_actions = input("Are you looking for a calculation to find the value of the coefficient of convection? If "
                            "so, then please tell us what method you want to use - (natural convection, or "
                            "forced convection). If not, then you must be here to calculate the heat transfer"
                            "in an enclosed space with an air gap. SO type anything else. : ")
    # creating an empty list once more
    list_of_dimensions = []
    # so, if they need to find the value of the coefficient of natural convection
    if further_actions == "natural convection":
        # we will first create a for-loop
        for i in range(4):
            # and ask them to put in the necessary data to find what they need exactly
            dimensions_needed = input("Please tell us the dimensions of the surface that will be affected by "
                                      "convection: (L, W, Diameter, and Theta: ")
            # appending the values into an empty list
            list_of_dimensions.append(float(dimensions_needed))
        # asking the user for the temperature of the air in case
        Tconvec = float(input("What is the temperature of the air/fluid?:  "))
        # now asking for the temperature of the surface just in case
        Ts = float(input("What is the temperature of the surface of the object? (Average it out if needed): "))

        # CALLING THE FUNCTION TO NOW FIND THE VALUE OF THE COEFFICIENT OF CONVECTION
        coeff_of_convection = natural_conv_h(list_of_dimensions[0], list_of_dimensions[1], list_of_dimensions[2],
                                             list_of_dimensions[3], values_returned[2], 9.81, 1/tfilm, values_returned[6],
                                             values_returned[5], Ts, Tconvec)
        # PRINTING OUT AND REVEALING THE CALCULATED VALUES OF h, Ra, AND Nu
        print(coeff_of_convection[0], coeff_of_convection[1], coeff_of_convection[2])
    # the other action possible is for enclosed surfaces
    else:
        # same process as shown above
        for i in range(3):
            # asking the user to now input the most vital dimensions of the window
            dimensions_needed = input("Please tell us the dimensions of the surface that will be affected by "
                                      "convection: (L, H, W: ")
            # creating a list from the values input
            list_of_dimensions.append(float(dimensions_needed))
        # reminding our user of how the program works
        print("So, we will once again use the lower and upper bounds needed so that this is AS accurate "
              "as possible.")
        # now, CALLING THE FUNCTION THAT WILL NOW GIVE US THE h VALUE OF THE ENCLOSED SURFACE FOR ACCURACY!!!
        coeff_of_convection = enclosure_checker(list_of_dimensions[0], list_of_dimensions[1], list_of_dimensions[2],
                                                values_returned[2], 9.81, values_returned[6], values_returned[5],
                                                1/tfilm, temp_1, temp_2)
        # printing out the values of h, Ra, and Nu
        print(coeff_of_convection[0], coeff_of_convection[1], coeff_of_convection[2])





