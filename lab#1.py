from math import sqrt, exp, log, ceil, floor
from sympy import primerange, legendre_symbol, factorint
import numpy as np

class rho_method():

    def __init__(self, n):
        self.n = n

    def rho_Pollard_factorization(self):
        x, y = 2, 2
        while True:
            x_new = self.calculate_xi(x, self.n)
            y_new = self.calculate_yi(y, self.n)
            d = self.calculate_gcd(x_new, y_new, self.n)
            x = x_new
            y = y_new
            if x == y:
                print('Алгоритм не знайшов дільника із заданими початковими значеннями')
                return 0
            if d != 1 and d != self.n:
                return d


    def calculate_f(self, x, n):
        return (x**2 + 1) % n

    def calculate_xi(self, x, n):
        return self.calculate_f(x, n)

    def calculate_yi(self, y, n):
        return self.calculate_f(self.calculate_f(y, n), n)

    def calculate_gcd(self, x, y, n):
        dif = x - y
        d = self.extended_EA(dif, n)[0]
        return d


    def extended_EA(self, a, b):
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


class BM_algorithm():

    def __init__(self, n):
        self.n = n
        self.a = 1/sqrt(2)
        self.factor_base = self.generate_factor_base()
        self.B_numbers = self.calculate_B_numbers()
        self.vectors = self.create_vectors()

    def generate_factor_base(self):
        size = ceil((exp(sqrt(log(self.n) * log(log(self.n))))) ** self.a)
        test_factor_base = [prime for prime in primerange(2, size+1)]
        final_factor_base = []
        for i in test_factor_base:
            if i == 2:
                final_factor_base.append(i)
                continue
            if legendre_symbol(n, i) == 1:
                final_factor_base.append(i)
        return final_factor_base
    
    def calculate_B_numbers(self):
        b_numbers = []
        alpha = sqrt(self.n)
        a = floor(alpha)
        v = 1
        u = a
        b_1 = 1
        b_2 = 0
        counter = 0
        while len(b_numbers) <= len(self.factor_base):
            v_new = int((n - u**2)/v)
            alpha_new = (sqrt(n) + u)/v_new
            a_new = floor(alpha_new)
            u_new = a_new * v_new - u
            b_new = a * b_1 + b_2
            b_new_sq = (b_new ** 2) % self.n
            can_form_b_new_sq = [i for i in factorint(b_new_sq).keys()]
            check = []
            for i in can_form_b_new_sq:
                if i in self.factor_base:
                    check.append(i)
            if check == can_form_b_new_sq:
                b_numbers.append(b_new_sq)
            v = v_new
            alpha = alpha_new
            a = a_new
            u = u_new
            b_1, b_2 = b_new, b_1
            counter += 1
        return b_numbers
    
    def create_vectors(self):
        vectors = []
        for i in range(len(self.B_numbers)):
            v = [0] * len(self.factor_base)
            b_num_main = [j for j in factorint(self.B_numbers[i]).keys()]
            b_num_power = [j for j in factorint(self.B_numbers[i]).values()]
            for k in range(len(b_num_main)):
                ind = self.factor_base.index(b_num_main[k])
                v[ind] = b_num_power[k] % 2
            vectors.append(v)
        return vectors
    
    def solve_system(self):
        vectors = self.vectors
        print(vectors)
        matrix_A = np.array(vectors).T.tolist()
        pass
        
n = 1495056764861639599
d = rho_method(n).rho_Pollard_factorization()
print(d)
n = int(n/d)

while rho_method(n).rho_Pollard_factorization() != 0:
    
    d = rho_method(n).rho_Pollard_factorization()
    n = int(n / d)
    print(f'n = {n}')
    print(f'd = {d}')