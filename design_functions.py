import thermodynamics as td
import massbalance as mb
import numpy as np

tau = np.array([
    [0, -0.18833332],
    [4.10588728, 0]
])
Psat_Water_bottom = td.Antoine_eqn(117.6, [3.55959	, 643.748, -198.043], True)
Psat_Butanol_bottom = td.Antoine_eqn(117.6, [4.50393,1313.878,-98.789], True)
Psat_Water_top = td.Antoine_eqn(93, [5.08354,1663.125,-45.622], True)
Psat_Butanol_top = td.Antoine_eqn(93, [4.50393,1313.878,-98.789], True)
def Fenske_multicomponent(T_D, T_B, x_D, x_B, Psat_i_dist, Psat_i_bot, Psat_j_dist, Psat_j_bot, tau):
    #x_D = (a, b, c, d..) mol basis = x_B
    gamma_D, _d = td.NRTL(x_D, T_D, tau)
    gamma_B, _b = td.NRTL(x_B, T_B, tau)
    alphaLKHK_D = (gamma_D[1] * Psat_i_dist) / (gamma_D[0] * Psat_j_dist)
    alphaLKHK_B = (gamma_B[1] * Psat_i_bot) / (gamma_B[0] * Psat_j_bot)
    geo_alpha = np.sqrt(alphaLKHK_B*alphaLKHK_D)
    print(alphaLKHK_B, Psat_i_bot, Psat_j_bot)
    print(alphaLKHK_D, Psat_i_dist, Psat_j_dist)
    Nmin = np.log((x_D[1]*x_B[0])/(x_D[0]*x_B[1])) / np.log(geo_alpha)
    print(Nmin)
    return

Fenske_multicomponent(93, 117.6, np.array([[0.2498],[0.7502]]), np.array([[0.3469], [0.653]]), Psat_Water_top, Psat_Water_bottom, Psat_Butanol_top, Psat_Butanol_bottom, tau)
