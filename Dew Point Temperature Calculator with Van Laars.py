"""
File:    Dew Point Temperature Calculator with Van Laars
Author:  Asad Siddiqui
Date:    11/6/2021
Section: Online Class, 11 AM
E-mail:  asiddiq3@umbc.edu
Description: The goal of this program is to figure out the best and somewhat optimal K values through a large scale
                iterative approach for the Dew Point Temperature, assuming the vapor molar composition is 50/50.
"""
import random


def recursive_kfactor(yi):
    # K1 and K2 will be randomized each time, and my rationality behind any where form 0 to 3 for both of them is that
    # the DePriester chart I found online has a low range for K2 which is the n-Pentane actually.
    K1 = random.uniform(0, 3)
    K2 = random.uniform(0, 3)

    # essentially, we will calculate x1 and x2 everytime using the general equation for K1 and K2. (just to make sense
    # of the values obtained and test for how much error there is within it).
    x1 = yi / K1
    x2 = yi / K2
    # we will sum up our x1 and x2 to see how much it is in total
    xtotal = x1 + x2
    # subtract from 1, and multiply by 100 to see how many percentage points it is off by
    error = (1 - xtotal) * 100
    # if the total liquid composition for some reason is greater than 1, that means it is not accurate AT ALL. And we
    # want a narrow range of liquid mole fractions with low error such that we can get an accurate iteration via
    # recursion
    if xtotal <= 1 and xtotal >= 0.97 and error < 3:
        # returning the K-factor values, and the respective liquid molar fractions for n-butane and n-pentane respective
        return (K1, K2, x1, x2)
    # otherwise,
    else:
        # print out the guessing of the K-factor values
        print(K1, K2)
        # and keep recursively finding the best values
        return recursive_kfactor(yi)


# main function
if __name__ == "__main__":
    # we were told the Z value was around 0.5, and so I assumed that if we were to find the Dew Point data, y1y2 are 0.5
    Y1 = 0.5
    Y2 = 0.5
    # the given pressure of the system
    P = 200  # (psi)
    # since both are the same, let's use one variable to describe the vapor mole fractions of the problem
    yi = Y1
    # print out the result from the recursive method
    print(recursive_kfactor(yi))
