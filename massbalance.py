import numpy as np
from scipy.optimize import root

R = 0.08206 #atm/J/k
Streams = {

}

def set_stream(Sn=str, S = 0, B = 0, A = 0, E = 0, W = 0, CO2 = 0, H2 = 0, NH4OH = 0,NH4_salts = 0):
    masses = [round(i, 2) for i in [S,B,A,E,W,CO2,H2,NH4OH,NH4_salts]]
    Streams[Sn] = masses
    return

def system_input(target_butanol):
    Temperature = 37 + 273.15
    sugar = target_butanol/0.18 #Y_btOH g/g
    water_sugar = sugar/0.06
    NH4OH = 35.046 * ( ( (sugar * 0.2 - target_butanol * 1/6) / 59.04) + ( (sugar * 0.25 - target_butanol) / 87.10) )
    water_base = 0.7 * NH4OH / 0.3
    v_broth = 72 * (((sugar + water_sugar) / 1.021) + ((NH4OH + water_base) / 0.9))
    CO2_volumetric = (sugar * 0.6 * R * Temperature) / 0.04401  # L/h
    CO2_stripping = (v_broth - CO2_volumetric) * (0.04401 / (R * Temperature) )  # kg/h
    set_stream("S1", S=sugar, W=water_sugar)
    set_stream("S2", W=water_base, NH4OH=NH4OH)
    set_stream("S19", CO2=CO2_stripping)
    R1(v_broth)
    return

def R1(v_broth):
    S1 = Streams["S1"]
    S2 = Streams["S2"]
    S19 = Streams["S19"]
    Total = sum(S1) + sum(S2) + sum(S19)
    # Totals
    BtOH_total = S1[0] * 0.3 * 0.6
    Actn_total = (BtOH_total * 10/6) * 0.3
    EtOH_total = (BtOH_total * 10/6) * 0.1
    CO2_fermentation = S1[0] * 0.6
    H2 = S1[0] * 0.0224
    CO2_dissolved = 0.033 * (v_broth/72) * 0.04401  # kg/h

    # Gas Stripping
    BAE_broth = [BtOH_total, Actn_total, EtOH_total]
    BAE_selectivities = [33.2, 12.1, 6.3]
    BAE_condensate_mass_fractions = []
    for i in range(len(BAE_broth)):
        x = BAE_broth[i]/Total
        y = (BAE_selectivities[i] * (x/(1-x)))/(1 + BAE_selectivities[i] * (x/(1-x)))
        BAE_condensate_mass_fractions.append(y)
    B = BtOH_total * 0.63
    A = B * (BAE_condensate_mass_fractions[1]/BAE_condensate_mass_fractions[0])
    E = B * (BAE_condensate_mass_fractions[2]/BAE_condensate_mass_fractions[0])
    W = (B + A + E) / sum(BAE_condensate_mass_fractions)
    set_stream("S3", B=B, A=A, E=E, W=W, H2=H2, CO2=(S19[5] + CO2_fermentation - CO2_dissolved))

    # Neutralization
    NH4_salts = ((( S1[0] * 0.2 - EtOH_total) / 59.04) * 77.083 + ((S1[0] * 0.25 - BtOH_total) / 87.10) * 105.137)
    water = S2[7] - NH4_salts

    # Waste Treatment
    BtOH_waste = BtOH_total - B
    Actn_waste = Actn_total - A
    EtOH_waste = EtOH_total - E
    sugar = S1[0] - sum([BtOH_total, Actn_total, EtOH_total, CO2_fermentation, H2])
    water = water + S1[4] + S2[4] - W
    set_stream("Waste", S=sugar, B=BtOH_waste, A=Actn_waste, E=EtOH_waste, W=water, CO2=CO2_dissolved, NH4_salts=NH4_salts)
    F1()

def F1():
    #
    S3 = Streams["S3"]
    set_stream("S4", S=S3[0], B=S3[1], A=S3[2], E=S3[3], W=S3[4], NH4OH=S3[7], NH4_salts=S3[8])
    set_stream("S18", CO2=S3[5], H2=S3[6])
    return

# Generic NRTL function
def NRTL(x, T, tau=None, a=None, b=None):
    N = x.shape[0] # Rows(Components)
    P = x.shape[1] # Columns(Phases)
    # alpha, tau, G = NxN and x, gamma = NxP
    alpha = np.full((N, N), 0.3)
    if tau is None:
        tau = a + b/T
        G = np.exp(-alpha * tau)
    elif a is None and b is None:
        G = np.exp(-alpha * tau)
    else:
        print("wrong input:\n if tau = None then a, b must be input\n if a,b = None then tau must be input")
    gamma = np.zeros((N, P))
    for phase in range(x.shape[1]):
        xi = x[:, phase]  # mol frac is the current phase
        term1 = np.zeros(N)
        term2 = np.zeros(N)
        for i in range(N):  # for each component(row) in the array tau (0, 1)
            term1[i] = np.sum(xi * tau[:, i] * G[:, i]) / np.sum(
                xi * G[:, i])  # term1 vector position is calculated by the column(j) of the current component(i)
            sum_j = 0
            for j in range(N):  # for every column repeated the number of times = items in xi (4)
                sum_mj = np.sum(xi * tau[:, j] * G[:,
                                                j])  # add up every x in the current phase * the columns of tau and G excluding the current one
                sum_kj = np.sum(xi * G[:, j])
                sum_j += (xi[j] * G[i, j] / sum_kj) * (tau[i, j] - sum_mj / sum_kj)
            term2[i] = sum_j
        ln_gamma = term1 + term2
        gamma[:, phase] = np.exp(ln_gamma)
    return gamma, tau

def composition_calculator(stream, temperature):
    T = temperature
    S = stream
    # Find Gamma BtOH/Water
    def unknown_components(T):
        # experimental LLE data
        xi = np.array([
            [0.488, 0.0191], #BtOH [org, aq]
            [0.512, 0.9809], #Water [org, aq]
        ])
        def nonlinear_regression(tau_vars):
            # initial guess
            tau = np.array([
                [0, tau_vars[0]], # Butanol(1) - Water(2)
                [tau_vars[1], 0]  # Water(2) - Butanol(1)
            ])
            gamma, tau = NRTL(xi, T, tau=tau)
            return gamma, tau
        def error_function(tau_vars):
            gamma, _ = nonlinear_regression(tau_vars)
            error = xi[:, 0] * gamma[:, 0] - xi[:, 1] * gamma[:, 1]
            return error
        #initial guess
        tau_guess = [1, 1]
        sol = root(error_function, tau_guess, method="hybr", tol=1e-3)
        if sol.success:
            optimal_tau_var = sol.x
            final_gamma, final_tau = nonlinear_regression(optimal_tau_var)
            error = nonlinear_regression(optimal_tau_var)
            return final_gamma, final_tau # BtOH, Water [org, aq]
        else:
            return None
    # Find Gamma EtOH/Actn/Water
    def known_components(S, T):
        # Relative Molar Fractions
        n1 = S[2] / 58.08
        n2 = S[3] / 46.069
        n3 = S[4] / 18.015
        nsum = n1 + n2 + n3
        xi = np.array([
            [n1/nsum],
            [n2/nsum],
            [n3/nsum]
        ])
        # Binary Interaction Parameters
        a = np.array([
            [0, 1.079, 6.398], #Acetone - 0
            [0.347, 0, 3.458], #Ethanol - 1
            [0.054, 0.801, 0]  #Water - 2
        ])
        b = np.array([
            [0, 479.1, 1809],  #Acetone - 0
            [206.6, 0, 586.1], #Ethanol - 1
            [420, 246.2, 0]    #Water - 2
        ])
        gamma, tau = NRTL(xi, T, a=a, b=b)
        return gamma, tau
    # Use Gamma to find xi' and xi"
    def LLE_solver():
        gamma_uc, tau_uc = unknown_components(T)
        gamma_kc, tau_kc = known_components(S, T)
        # Interaction parameter tau matrix
        BA = 331.72/T
        BE = 173.2/T
        AB = -120.47/T
        EB = -90.84/T
        tau = np.array([
            [0, BA, BE, tau_uc[0,1]],  # B
            [AB, 0, tau_kc[0,1], tau_kc[0,2]],  # A
            [EB, tau_kc[1,0], 0, tau_kc[1,2]],  # E
            [tau_uc[1,0], tau_kc[2,0], tau_kc[2,1], 0]  # W
        ])
        # Stream Composition zi
        nB = S[1] / 74.123
        nA = S[2] / 58.08
        nE = S[3] / 46.069
        nW = S[4] / 18.015
        nsum = nB + nA + nE + nW
        zi = np.array([
            [nB/nsum], [nA/nsum], [nE/nsum], [nW/nsum]
        ])
        N = zi.shape[0]
        phases = 2
        for P in range(phases):
            pass
        return
    # Assume EtOH / BtOH, Actn / BtOH fully miscible therefore gamma = 1
    LLE_solver()
    return

def Dec1():
    S4 = Streams["S4"]
    Total = sum(S4)
    w_BtOH_in = S4[1]/Total
    w_Water_in = S4[4]/Total
    # Lever rule - Dec1 layer composition
    C_aq = 0.0742 # empirical - BtOH in aqueous
    C_org = 0.797 # empirical - BtOH in organic
    w_aq = (C_org - w_BtOH_in) / (C_org - C_aq)
    w_org = (w_BtOH_in - C_aq) / (C_org - C_aq)
    m_org = Total * w_org
    m_aq = Total * w_aq
    #Assumption : BtOH distribution only relies on these weight fractions, the other components use NRTL model and the salts all remain in the aqueous layer.
    #Organic Layer
    B_org = C_org * S4[1]

    #Aqueous Layer
    B_aq = C_aq * S4[1]

system_input(63.13)
composition_calculator(Streams['S4'], 310)

"""
print("kg/h   S, B, A, E, W, CO2, H2, NH4OH, Salts")
for a in Streams:
    print(f"{a}, {Streams[a]} = {sum(Streams[a])}")

print(f"mass balance: in-out {sum(Streams['S3'])  } - {sum(Streams['S4']) + sum(Streams['S18'])}")"""