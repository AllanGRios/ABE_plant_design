import matplotlib.pyplot as plt
import thermodynamics as td
import numpy as np

def Txy_draw(P, T_lower, T_upper, components):
    # Pressure in kPa, Temp in Kelvin
    # Components BtOH, Actn, EtOH, Water. index 0 -> 3
    T, x, y = td.Txy(P, T_lower, T_upper, components)
    plt.plot(x, T, label="X - Bubble Point")
    plt.plot(y, T, label="Y - Dew Point")
    plt.xlabel(f"X {components[0]}")
    plt.ylabel("Temperature")
    plt.show()
    return

# Raoult law doesn't work for non ideal mixtures