import sympy as sym

sym.init_printing()

Z, P, V, T, R, Tc, Pc, a, b, A, B = sym.symbols('Z, P, V, T, R, Tc, Pc, a, b, A, B')

P = 1500000
T = 323.15
R = 8.314
Tc = 190.56
Pc = 4599000

a = 0.42798 * (((R ** 2) * (Tc ** 2.5)) / Pc)
b = 0.08664 * (R * Tc / Pc)

A = (a * P) / ((R ** 2) * (T ** 1.5))
B = (b * P) / ((R) * (T))

volume = sym.Eq(P - (((R * T) / (V - b)) -
                               (a / ((T ** 0.5) * V * (V + b)))))

compressibility_factor_Z = sym.Eq((Z ** 3) - (Z ** 2) + (A - B - (B ** 2)) * Z - (A * B))

print(sym.solve(volume, V))

print(sym.solve(compressibility_factor_Z, Z))