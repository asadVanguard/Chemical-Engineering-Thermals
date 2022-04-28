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
    critical_temperature = float(input('What is the value of the critical temperature? '))

    critical_pressure = float(input('What is the value of the critical pressure? '))

    func_a = 0.42798 * (((constant_r ** 2) * (critical_temperature ** 2.5)) / critical_pressure)

    func_b = 0.08664 * (constant_r * critical_temperature / critical_pressure)



