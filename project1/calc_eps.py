from decimal import Decimal as dec

def eps(n):
    h = dec(1) / dec(1+n)
    H = dec(10)*h
    return (dec(1) + dec(100)*dec(h)*dec(h)/(2 - H.exp() - (-H).exp())).log10()

if __name__ == "__main__":
    for i in range(8):
        n = 10**i
        print(f"{n}, {eps(n):.7f}")
