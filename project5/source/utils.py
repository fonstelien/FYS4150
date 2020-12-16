import numpy as np
import sklearn as skl
from sirs import Sirs, SirsSolver

def make_sirs_list(n, *args, **kwargs):
    sirs_list = list()
    for i in range(n):
        sirs_list.append(Sirs(*args, **kwargs))
    return sirs_list

def mc_stats(sirs_list, t=20, seed=0):
    n = len(sirs_list)
    S = np.ndarray(n)
    I = np.ndarray(n)
    R = np.ndarray(n)
    N = np.ndarray(n)
    
    for i, sirs in enumerate(sirs_list):
        solver = SirsSolver(sirs)
        solver.run_mc(t_max=t, seed=i+seed)
        S[i] = sirs.S
        I[i] = sirs.I
        R[i] = sirs.R
        N[i] = sirs.N

    res = {'S':(S.mean(), S.std()),
           'I':(I.mean(), I.std()),
           'R':(R.mean(), R.std()),
           'N':(N.mean(), N.std())}

    return res
        
