import numpy as np
from scipy.optimize import root
from scipy.optimize import least_squares

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

def LLE_solver(zi, tau, T, xi_guess, beta_guess, MW=None):
    N = zi.shape[0]
    unknowns = N*2 + 1
    #initial flattening
    xi_guess = xi_guess.flatten()
    x0_guess = np.concatenate([xi_guess, [beta_guess]])

    # Minimize error to find root
    def error_function_LLE(x0_guess):
        xi = x0_guess[:N*2].reshape(N, 2)
        beta = x0_guess[N*2]
        gamma, _ = NRTL(xi, T, tau=tau)
        isoactivity_error = xi[:, 0] * gamma[:, 0] - xi[:, 1] * gamma[:, 1]  # x0' * γ' - x0" * γ" = ε
        mol_balance_error = zi[:, 0] - beta * xi[:, 0] - (1 - beta) * xi[:, 1]
        mass_balance_error_org = np.array([np.sum(xi[:, 0]) - 1])
        error = np.concatenate([isoactivity_error, mol_balance_error, mass_balance_error_org])
        return error

   # least squares method solver
    lower = np.zeros(unknowns)
    upper = np.concatenate([np.ones(N*2), [1]])
    sol = least_squares(error_function_LLE, x0_guess, bounds=(lower, upper))
    if sol.success:
        final_xi = sol.x[:8].reshape(N, 2)
        final_beta = sol.x[8]
        if MW is not None:
            mass_org = final_xi[:, 0] * MW
            w_org = mass_org / np.sum(mass_org)
            mass_aq = final_xi[:, 1] * MW
            w_aq = mass_aq / np.sum(mass_aq)
            return final_xi, final_beta, w_org, w_aq
        else:
            final_gamma, final_tau = NRTL(final_xi, T, tau=tau)
            return final_xi, final_beta, final_gamma
    else:
        return
    return

def fenske_multicomponent(zi=None, tau=None, T=298, yx_guess=None, beta_guess=None, MW=None):
    #if tau is not None:
    tau = np.array([
        [0, 2, 3, 4],
        [5, 0, 7, 3],
        [9, 10, 0, 2],
        [2, 5, 6, 0]
    ])
    zi = np.array([
        [2],
        [3],
        [8],
        [9]
    ])
    yx_guess = np.array([
        [0.3, 0.02],
        [0.2, 0.6],
        [0.4, 0.08],
        [0.1, 0.2]
    ])
    beta_guess = 0.07
    N = zi.shape[0] # Components in feed
    P = 2 # Phases (Vapor Liquid)
    """
     1. VLE NRTL get gamma values
     2. Calculate relative volatitlity value for LK HK
     3. calculate N+1 stages using fenske assuming total reflux"""
    yx,  beta, final_gamma = LLE_solver(zi, tau, T, yx_guess, beta_guess)
    print(yx)
    print(beta)
    print(final_gamma)
    #alpha = (gamma_LK * Pvap_LK) / (gamma_HK * Pvap_HK)
    return

fenske_multicomponent()