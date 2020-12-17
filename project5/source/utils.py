import numpy as np
import pandas as pd

from sirs import Sirs, SirsSolver


def exponential_decay(t, lmd, floor=0., delay=0.):
    '''Exponential decay function. Keeps value 1 for t<delay; then exponential decay with time constant (half life time) lmd lmd; converges to floor for t>delay.'''
    if t < delay:
        return 1
    tau = lmd/np.log(2)
    return (1-floor)*np.exp(-1*(t-delay)/tau) + floor


def exponential_growth(t, lmd, floor=0., delay=0.):
    '''Exponential growth function. Grows exponentially from floor until value 1 at t=delay with time constant (half life time) lmd lmd; then keeps value 1.'''
    if t > delay:
        return 1
    tau = lmd/np.log(2)
    return (1-floor)*np.exp((t-delay)/tau) + floor


def bathtub(t, lmd, floor=0, delay=0, duration=10):
    '''Combination of exponential_decay(t, lmd, delay=dealy, floor=floor) and exponential_growth(t, lmd, floor=floor, delay=delay+duration).'''
    decay = exponential_decay(t, lmd, floor=floor, delay=delay)
    growth = exponential_growth(t, lmd, floor=floor, delay=delay+duration)
    return decay + growth - floor


def seasonal_variation(t, a=.5, A=.5, multiplier=1, freq_multiplier=1, shift=0):
    '''Seasonal variation with mean a*multiplier, amplitude A*multiplier, period 10*freq_multiplier and phase shifted by shift.'''
    a, A = a*multiplier, A*multiplier
    freq = 2*np.pi/10 * freq_multiplier
    return a + A*np.cos(freq*(t - shift))


def make_sirs_list(n, *args, **kwargs):
    '''Returns a list of n identical Sirs instances instatiated as Sirs(*args, **kwargs).'''
    sirs_list = list()
    for i in range(n):
        sirs_list.append(Sirs(*args, **kwargs))
    return sirs_list


def mc_stats(sirs_list, seed, *args, **kwargs):
    '''Runs Monte Carlo simulation on each Sirs instance in sirs_list. Calculates the mean and standard deviations of each Sirs population Sirs.S,I,R,N at over time and returns them as time series (mean,std) in a dict() keyed on "S","I","R","N".'''
    n = len(sirs_list)
    S_df = pd.DataFrame(index=[.0])
    I_df = S_df.copy()
    R_df = S_df.copy()
    N_df = S_df.copy()
    
    for i, sirs in enumerate(sirs_list):
        solver = SirsSolver(sirs)
        solver.run_mc(*args, seed=seed+i, **kwargs)
        S, I, R, N = solver.get_fractions()
        S_df = S_df.join(pd.DataFrame(data=S, index=solver.t, columns=[i]), how='outer')
        I_df = I_df.join(pd.DataFrame(data=I, index=solver.t, columns=[i]), how='outer')
        R_df = R_df.join(pd.DataFrame(data=R, index=solver.t, columns=[i]), how='outer')
        N_df = N_df.join(pd.DataFrame(data=N, index=solver.t, columns=[i]), how='outer')

    S_df = S_df.ffill().bfill()
    I_df = I_df.ffill().bfill()
    R_df = R_df.ffill().bfill()
    N_df = N_df.ffill().bfill()
    S_mean = S_df.mean(axis=1).to_numpy()
    I_mean = I_df.mean(axis=1).to_numpy()
    R_mean = R_df.mean(axis=1).to_numpy()
    N_mean = N_df.mean(axis=1).to_numpy()
    S_std = np.sqrt((n/(n-1))*S_df.var(axis=1)).to_numpy()
    I_std = np.sqrt((n/(n-1))*I_df.var(axis=1)).to_numpy()
    R_std = np.sqrt((n/(n-1))*R_df.var(axis=1)).to_numpy()
    N_std = np.sqrt((n/(n-1))*N_df.var(axis=1)).to_numpy()

    t = S_df.index
    res = {'S':(S_mean, S_std),
           'I':(I_mean, I_std),
           'R':(R_mean, R_std),
           'N':(N_mean, N_std)}

    return t, res
