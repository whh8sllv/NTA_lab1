from math import sqrt, exp, log, ceil, floor
from sympy import primerange, legendre_symbol, factorint


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
                return 'Алгоритм не знайшов дільника із заданими початковими значеннями'
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
    
    def B_numbers(self):
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
            print('*********')
            print(f'v = {v_new}')
            print(f'alpha = {alpha_new}')
            print(f'a = {a_new}')
            print(f'u = {u_new}')
            print(f'b = {b_new}')
            print(f'b_new_sq = {b_new_sq}')
            print('*********')
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
        return self.factor_base, b_numbers
    
    

            

            


    



n = 17873

div1 = rho_method(n).rho_Pollard_factorization()
print(div1)

div2 = BM_algorithm(n).B_numbers()
print(div2)