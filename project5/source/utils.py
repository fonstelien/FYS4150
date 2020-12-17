import numpy as np
import sklearn as skl
from threading import Thread

from sirs import Sirs, SirsSolver


def make_sirs_list(n, *args, **kwargs):
    '''Returns a list of n identical Sirs instances instatiated as Sirs(*args, **kwargs).'''
    sirs_list = list()
    for i in range(n):
        sirs_list.append(Sirs(*args, **kwargs))
    return sirs_list


def mc_stats(sirs_list, t=20, seed=0):
    '''Runs Monte Carlo simulation up to time t on each Sirs instance in sirs_list. Calculates the mean and standard deviations of each Sirs population Sirs.S,I,R,N at time t and returns them in a dict() keyed on "S","I","R","N" along with the samples; (S,I,R,N).'''
    n = len(sirs_list)
    S = np.ndarray(n)
    I = np.ndarray(n)
    R = np.ndarray(n)
    N = np.ndarray(n)

    for i, sirs in enumerate(sirs_list):
        solver = SirsSolver(sirs)
        solver.run_mc(t_max=t, seed=seed+i)
        S[i] = sirs.S
        I[i] = sirs.I
        R[i] = sirs.R
        N[i] = sirs.N
        
    res = {'S':(S.mean(), S.std()),
           'I':(I.mean(), I.std()),
           'R':(R.mean(), R.std()),
           'N':(N.mean(), N.std())}
    
    return res, (S, I, R, N)
