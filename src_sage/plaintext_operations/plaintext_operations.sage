Z = PolynomialRing(ZZ, 'x') # integers R1
Q = PolynomialRing(QQ, 'x') # rationals K1

K.<x, y, z> = PolynomialRing(QQ) # rationals in K

from sage.structure.element import is_Vector        
class FieldK:
    def __init__(self, prime_powers):
        assert(len(prime_powers) == 3)
        assert(prime_powers[0][0] != prime_powers[1][0] != prime_powers[2][0])
        for i in range(3):
             assert(is_prime(prime_powers[i][0]))
        
        self.prime_powers = prime_powers
        self.pm = [prime_powers[i][0] ** prime_powers[i][1] for i in range(3)]
        self.pm.append(self.pm[0] * self.pm[1] * self.pm[2])
        self.r = min(euler_phi(self.pm[1]), euler_phi(self.pm[2]))
        self.cyc = Z(cyclotomic_polynomial(self.pm[0]))
        self.ncyc = Z(cyclotomic_polynomial(self.pm[3]))

    def random_element(self, mode = -1):
        if mode == -1:
            mode = ZZ.random_element(1, 4)
        new_polys = vector([Z.random_element(degree = euler_phi(self.pm[0]) - 1) for _ in range(self.r)])
        
        return Plaintext(new_polys, mode, self.cyc, self.r)
        
    def create_element(self, polys, mode):
        return Plaintext(polys, mode, self.cyc, self.r)

    def compose_poly(self, plaintext):
        result = Z(0)
        assert(plaintext.mode == 1 or plaintext.mode == 2)
        if plaintext.mode == 1:
            for i in range(self.r):
                result = result + (Z(x^(i * self.pm[0])) * plaintext.polys[i](x = x^(self.pm[1]))) % cyclotomic_polynomial(self.pm[0] * self.pm[1]) # R12

        if plaintext.mode == 2:
            for i in range(self.r): 
                result = result + (Z(x^(i * self.pm[0])) * plaintext.polys[i](x = x^(self.pm[2]))) % cyclotomic_polynomial(self.pm[0] * self.pm[2]) # R13

        return result
        
class Plaintext: 
    def __init__(self, polys, mode, cyc, r):
        assert(0 < mode < 5)
        assert(is_Vector(polys))
        assert(len(polys) == r) 

        for j in range(len(polys)):
            polys[j] = Z(polys[j]) % cyc

        self.mode = mode
        self.polys = polys
        self.cyc = cyc 
        self.r = r


    def __add__(self, other):
        assert(self.mode == other.mode)
        new_polys =  self.polys + other.polys
        return Plaintext(new_polys, self.mode, self.cyc, self.r)

    def __mul__(self, other):   
        assert(abs(self.mode - other.mode) == 2)

        if (self.mode, other.mode) in [(1,3),(3,1)]:
            mode = 2
        elif (self.mode, other.mode) in [(2,4),(4,2)]:
            mode = 1
        
        new_polys = self.polys.pairwise_product(other.polys)
        return Plaintext(new_polys, mode, self.cyc, self.r)

    def print_poly(self):
        assert(self.mode == 1 or self.mode == 2)

        result = 0
        if self.mode == 1:
            for i in range(self.r):
                result = result + K(y^i) * K(self.polys[i])
        elif self.mode == 2:
            for i in range(self.r):
                result = result + K(z^i) * K(self.polys[i])
        
        print(result)

'''
# Testing

# Define prime powers
p1, m1 = 2, 2
p2, m2 = 3, 2
p3, m3 = 5, 1

prime_powers = [(p1, m1),(p2, m2),(p3, m3)]

K_ = FieldK(prime_powers)

epochs = 10
for i in range(epochs):
    for pair in [(1,3),(2,4)]:
        pol1 = K_.random_element(pair[0])
        pol2 = K_.random_element(pair[1])

        print("polinomio 1 ",pol1.polys)
        print("polinomio 2 ",pol2.polys)
        pol1.print_poly()
        print("pol1 x pol2 ",(pol1 * pol2).polys)
'''