import matplotlib as plt
import thermodynamics as td
import numpy as np

tau = np.array([
    [0, 0.5],
    [0.7, 0]
])
P_sat = 893 #pa
def Txy_draw(P, P_sat, tau):
    temperature = np.linspace(273, 473)
    mol_frac = np.array([
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1
    ])
    for T in temperature:
        gamma,_ = NRTL(mol_frac, T, tau=tau)
