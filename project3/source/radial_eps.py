from decimal import Decimal as dec
from decimal import getcontext
import numpy as np

def eps_v(x0, h):
    x0 = dec(x0)
    h = dec(h)
    pi = dec(np.pi)
    v0 = 2*pi*(2*(1/x0-1)).sqrt()
    x1 = x0 + h*v0 + 2*(h*pi/x0)**2
    v1_true = 2*pi*(2*(1/x1-1)).sqrt()
    v1_verlet = v0 - 2*h*pi**2*(1/x1**2+1/x0**2)
    return abs((v1_true-v1_verlet)/v1_true), v1_true, v1_verlet, v0

if __name__ == "__main__":
    getcontext().prec = 50
    h = 1.0E-6
    n = 1000
    x = np.linspace(.001, 1-10*h, n)
    err = [eps_v(xi, h) for xi in x]
    print("x,eps,v_true,v_verlet,v0")
    for xi in x:
        eps, v_true, v_verlet, v0 = eps_v(xi, h)
        print(f"{xi},{eps},{v_true},{v_verlet},{v0}")
        
        
