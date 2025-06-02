def rho_Pollard_factorization(n):
    x, y = 2, 2
    d = 1
    for i in range(10):
        x_new = calculate_xi(x, n)
        y_new = calculate_yi(y, n)
        d = calculate_gcd(x_new, y_new, n)
        x = x_new
        y = y_new
        print(f'x new = {x}')
        print(f'y_new = {y}')


def calculate_f(x, n):
    return (x**2 + 1) % n

def calculate_xi(x, n):
    return calculate_f(x, n)

def calculate_yi(y, n):
    return calculate_f(calculate_f(y, n), n)

def calculate_gcd(x, y, n):
    dif = x - y
    d = extended_EA(dif, n)[0]
    return d


def extended_EA(a, b):
    u0, u1 = 1, 0
    v0, v1 = 0, 1
    r = a % b
    q = a // b
    while r != 0:
        u0, u1 = u1, (u0 - q * u1)
        v0, v1 = v1, (v0 - q * v1)
        a = b
        b = r
        r = a % b
        q = a // b
    return (b, u1, v1)


n = 527

rho_Pollard_factorization(n)