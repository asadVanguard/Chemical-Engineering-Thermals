import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
import sympy as sym
sym.init_printing()

if __name__ == '__main__':

    P, R, T, Z, Vm, a, b, A, B = sym.symbols('P, R, T, Z, Vm, a, b, A, B')

    Tcrit = 647.3
    Pcrit = 221.2
    T = 373.15
    Treduce = (T/Tcrit)

    w = 0.344
    R = 83.14

    b = (0.07780*R)*(Tcrit/Pcrit)
    k = 0.37469+(1.54226*w) - (0.2669*(w ** 2))

    alpha = (1 + (k*(1-(Treduce**0.5))))**2
    Ac = 0.45724 * ((R**2) * (Tcrit**2)/Pcrit)
    a = Ac*alpha


    sentinel = 1
    pressure_array = []
    while(sentinel != 0):
        pressures = float(input("What pressures do you want to find the fugacity at? (Press 0 to stop): "))
        if pressures == 0:
            sentinel = 0
        else:
            pressure_array.append(pressures)

    fugacity_array = []
    for i in range(len(pressure_array)):
        A = (a * pressure_array[i]) / (R ** 2 * T ** 2)
        B = (b * pressure_array[i]) / (R * T)

        PR_eqn = sym.Eq(pressure_array[i] - (((R * T) / (Vm - b)) - (a / ((Vm * (Vm + b)) + (b * (Vm - b))))))
        Volume = sym.solve(PR_eqn, Vm)
        PR_Z_eqn = sym.Eq(Z - (Volume[0] / (Volume[0] - b)) - ((a * Volume[0]) / (R * T)) * (
                    1 / ((Volume[0] * (Volume[0] + b)) + (Volume[0] * (Volume[0] - b)))))
        Z_value = sym.solve(PR_Z_eqn, Z)
        Z_converted = 1 - abs((Z_value[0] * (10 ** -2)))
        fug = pressure_array[i] * math.exp((Z_converted - 1) - (sym.log(Z_converted - B) - (A / (2 * B * (2 ** (1 / 2))))) * sym.log((Z_converted + (1 + 2 ** (1 / 2)) * B) / (Z_converted + (1 - 2 ** (1 / 2)) * B)))

        fugacity_array.append(fug)
        print(fugacity_array)
        print(fug)
    plt.plot(pressure_array, fugacity_array)
    plt.xlabel("Pressure (Bar)")
    plt.ylabel("Fugacities (Bar)")
    plt.title("Fugacity vs Pressure w/ Peng Robinson")
    plt.xlim(5, 150, 10)
    plt.ylim(5, 300, 10)
    plt.show()