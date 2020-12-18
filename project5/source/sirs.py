import numpy as np


class Sirs:
    def __init__(self, N0, I0, A=0, a=4, b=2, c=.5, d=0, dI=0, e=0, F=0, f=0, fa=None, ff=None):
        '''Instantiate with initial pop. size and number of infected I0. A,a are the amplitude and average of the seasonal infection rate; b is the recovery rate; c is the immunity loss rate; d is the population death rate from other causes; dI is the insidence death rate; e is the population birth rate; F,f are the amplitude and average of the seasonal vaccination rate; fa,ff accept functions func(t) for infection and vaccination rates.'''
        self.N = N0
        self.S = N0 - I0
        self.I = I0
        self.R = 0
        self.fa = fa if fa else (lambda t: A*np.cos(2*np.pi/10*t) + a)
        self.b = b
        self.c = c
        self.d = d
        self.dI = dI
        self.e = e
        self.ff = ff if ff else (lambda t: F*np.cos(2*np.pi/10*t) + f)
        
    def S_dot(self, ti, Si):
        '''S time derivative.'''
        return self.c*self.R - self.fa(ti)*Si*self.I/self.N - self.d*Si + self.e*self.N - self.ff(ti)*Si
    
    def I_dot(self, ti, Ii):
        '''I time derivative.'''
        return self.fa(ti)*self.S*Ii/self.N - (self.b + self.d + self.dI)*Ii
    
    def R_dot(self, ti, Ri):
        '''R time derivative.'''
        return self.b*self.I - (self.c + self.d)*Ri + self.ff(ti)*self.S
    
    def N_dot(self, ti, Ni):
        '''N time derivative.'''
        return self.e*Ni - self.d*self.N - self.dI*self.I
    
    def delta_t(self, ti=0):
        '''Returns time delta for Monte Carlo simulation.'''
        a = self.fa(ti)
        f = self.ff(ti)
        return min(4/(a*self.N) if a else 1, 
                   1/(self.b*self.N),
                   1/(self.c*self.N), 
                   1/(self.d*self.N) if self.d else 1, 
                   1/(self.dI*self.N) if self.dI else 1,
                   1/(self.e*self.N) if self.e else 1,
                   1/(f*self.N) if f else 1)
    
    def sir_transitions(self, random_number, dt, ti=0):
        '''Returns S->I->R->S transitions as (si, ir, rs, sr), where sr is the vaccination transitions.'''
        a = self.fa(ti)
        f = self.ff(ti)
        si = 1 if random_number < (a*self.S*self.I/self.N)*dt else 0  # S -> I
        ir = 1 if random_number < (self.b*self.I)*dt else 0           # I -> R
        rs = 1 if random_number < (self.c*self.R)*dt else 0           # R -> S
        sr = 1 if random_number < (f*self.S)*dt else 0                # S -> R by vaccination
        return (si, ir, rs, sr)
    
    def vitality_transitions(self, random_number, dt):
        '''Returns the transitions into and out of life as (bs, sd, idI, idD, rd), where bs are births; sd are deaths among the susceptible; idI are deaths from disease; idD are deaths from other causes among the infected; rd are deaths among the recovering.'''
        bs = 1 if random_number < (self.e*self.N)*dt else 0           # born -> S
        sd = 1 if random_number < (self.d*self.S)*dt else 0           # S -> dead
        idI = 1 if random_number < (self.dI*self.I)*dt else 0         # I -> dead disease
        idD = 1 if random_number < (self.d*self.I)*dt else 0          # I -> dead other cause
        rd = 1 if random_number < (self.d*self.R)*dt else 0           # s -> dead
        return (bs, sd, idI, idD, rd)
        
    def update(self, new_vals):
        '''Updates Sirs.S,I,R,N to new_values=(S_new, I_new, R_new, N_new).'''
        S, I, R, N = new_vals
        self.S = S
        self.I = I
        self.R = R
        self.N = N

    def get_fractions(self):
        N = self.N
        return (self.S/N, self.I/N, self.R/N)
        
    def constant_fractions(self):
        '''Returns constant fractions in a zero deaths/zero births population with constant transmission rate fa(t)=a.'''
        a = self.fa(0)
        f = self.ff(0)
        if a == 0 or a <= self.b:
            s = 1/(1 + f/self.c)
            r = 1 - s
            return (s,0,r)
        s = self.b/a
        i = (1 - (1+f/self.c)*self.b/a)/(1 + self.b/self.c)
        r = 1 - s - i
        return (s, i, r)


class SirsSolver:
    def __init__(self, sirs):
        '''Class for solving the SIRS model Instantiate with a Sirs instance.'''
        self.sirs = sirs
        self.fs = {'S':sirs.S_dot, 'I':sirs.I_dot, 'R':sirs.R_dot, 'N':sirs.N_dot}
        self.res = {'S':None, 'I':None, 'R':None, 'N':None}
        self.t = self.h = None

    def get_fractions(self):
        '''Returns the simulation results SirsSolver.res as fractions of dynamical population.'''
        N = self.res['N']
        return (self.res['S']/N, self.res['I']/N, self.res['R']/N, N)
        
    def _init_run(self, n_iter):
        '''Initializes the SirsSolver.res data structure. Do this before simulation.'''
        self.res['S'] = np.ndarray(n_iter)
        self.res['I'] = np.ndarray(n_iter)
        self.res['R'] = np.ndarray(n_iter)
        self.res['N'] = np.ndarray(n_iter)
        self.res['S'][0] = self.sirs.S
        self.res['I'][0] = self.sirs.I
        self.res['R'][0] = self.sirs.R
        self.res['N'][0] = self.sirs.N
    
    def _rk4(self, code, i):
        '''Implementation of RK4'''
        f = self.fs[code]
        yi = self.res[code][i-1]
        ti = self.t[i]
        h = self.h
        
        k1 = h*f(ti, yi)
        k2 = h*f(ti+.5*h, yi+.5*k1)
        k3 = h*f(ti+.5*h, yi+.5*k2)
        k4 = h*f(ti+h, yi+k3)
        
        return max(0, yi + (k1 + 2*k2 + 2*k3 + k4)/6)

    def run_rk4(self, t_max=20, n_iter=1000):
        '''Run 4th order Runge-Kutta simulation until t_max in n_iter iterations. Result is available in the SirsSolver.res data structure.'''
        self._init_run(n_iter+1)
        self.t = np.linspace(0, t_max, n_iter+1)
        self.h = t_max/n_iter
        for i in range(1, n_iter+1):
            self.res['S'][i] = self._rk4('S', i)
            self.res['I'][i] = self._rk4('I', i)
            self.res['N'][i] = self._rk4('N', i)
            self.res['R'][i] = self.res['N'][i] - self.res['S'][i] - self.res['I'][i]
            new_vals = (self.res['S'][i], self.res['I'][i], self.res['R'][i], self.res['N'][i])
            self.sirs.update(new_vals)
    

    def run_mc(self, t_max=20, seed=0):
        '''Run Monte Carlo simulation until t_max with seed forwarded to RNG. Step length is adaptive. Result is available in the SirsSolver.res data structure.'''
        n_iter = 10000
        self._init_run(n_iter)
        self.t = np.zeros(n_iter)
        np.random.seed(seed)
        random_numbers = np.random.rand(n_iter)
    
        ## Monte Carlo
        self.t[0] = ti = 0
        dt = self.sirs.delta_t(0)
        i = 1
        while ti <= t_max:
            if i%n_iter == 0:
                random_numbers = np.random.rand(n_iter)
                self.t.resize(i+n_iter, refcheck=False)
                self.res['S'].resize(i+n_iter, refcheck=False)
                self.res['I'].resize(i+n_iter, refcheck=False)
                self.res['R'].resize(i+n_iter, refcheck=False)
                self.res['N'].resize(i+n_iter, refcheck=False)
            
            random_number = random_numbers[i%n_iter]
            si, ir, rs, sr = self.sirs.sir_transitions(random_number, dt, ti=ti)
            bs, sd, idI, idD, rd = self.sirs.vitality_transitions(random_number, dt)
            
            self.t[i] = ti + dt
            self.res['S'][i] = self.res['S'][i-1] - si + rs + bs - sd - sr
            self.res['I'][i] = self.res['I'][i-1] + si - ir - idI - idD
            self.res['N'][i] = self.res['N'][i-1] + bs - sd - idI - idD - rd
            self.res['R'][i] = self.res['N'][i] - self.res['S'][i] - self.res['I'][i]
            new_vals = (self.res['S'][i], self.res['I'][i], self.res['R'][i], self.res['N'][i])
            self.sirs.update(new_vals)
            
            ti = self.t[i]
            dt = self.sirs.delta_t(ti)
            i += 1
            
        ## Trim trailing zeros
        self.t = np.trim_zeros(self.t, trim='b')
        new_length = len(self.t)
        self.res['S'] = self.res['S'][:new_length]
        self.res['I'] = self.res['I'][:new_length]
        self.res['R'] = self.res['R'][:new_length]
        self.res['N'] = self.res['N'][:new_length]


    def get_stats(self, t_start):
        '''Calculates the mean and standard deviations of each population S,I,R,N after time t_start and returns them in a dict() keyed on "S","I","R","N". To be used after run_mc().'''
        idx = self.t[self.t < t_start].size
        res = {'S':(self.res['S'][idx:].mean(), self.res['S'][idx:].std()),
               'I':(self.res['I'][idx:].mean(), self.res['I'][idx:].std()),
               'R':(self.res['R'][idx:].mean(), self.res['R'][idx:].std()),
               'N':(self.res['N'][idx:].mean(), self.res['N'][idx:].std())}
        return res, (self.res['S'][idx:], self.res['I'][idx:], self.res['R'][idx:], self.res['N'][idx:])
