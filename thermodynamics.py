import numpy as np

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