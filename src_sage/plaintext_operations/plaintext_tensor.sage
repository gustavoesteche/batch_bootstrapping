Q = PolynomialRing(QQ, 'x')
Z = PolynomialRing(ZZ, 'x')

load("./dual/trace_wdual")

#### Previously
# Calculating Trace R1 x R2 x R3 / R1 x R(2|3) 
def trace_mode(poly_m, factors_m, mode, index = 0):
    '''
    Computes Tr_K/K13 mode = 1 or Tr_K/K12 mode = 2
    index: which polynomial will be multiplied by the trace factor
    factors_m: A tuple list representing the prime power that divides m
    poly_m[i]: the polynomial in the field K_{p_i ^m_i} that decomposes the original poly 
    '''
    assert(len(poly_m) == 3), "3 fields are necessary"
    assert(mode == 1 or mode == 2), "Choose mode 1 or 2"

    trace = fast_trace_pm_to_1(factors_m[mode][0], factors_m[mode][1], poly_m[mode])
    poly_m[mode] = None
    poly_m[index] = poly_m[index] * trace 
    return poly_m
####

class TesorMField:
    def __init__(self, prime_powers):
        # assert there are 3 fields
        assert(len(prime_powers) == 3)
        
        # assert primes are different
        assert(prime_powers[0][0] != prime_powers[1][0] != prime_powers[2][0])
        
        # assert the numbers are really prime
        for i in range(3):
             assert(is_prime(prime_powers[i][0]))
        
        self.prime_powers = prime_powers
        self.cyclo_poly = [Z(cyclotomic_polynomial(prime_powers[i][0] ** prime_powers[i][1])) for i in range(3)] 
    
        
    def random_element(self, mode = -1):
        if mode == -1:
            mode = randint(1,4)  
        else:
            assert(0 < mode < 5)
            
        poly_m = [Z.random_element(degree=euler_phi(self.prime_powers[i][0] ** self.prime_powers[i][1])-1) for i in range(3)]
        
        if mode == 1:
            poly_m[2] = None
        elif mode == 2:
            poly_m[1] = None 
        elif mode == 3:
            poly_m[1] = sample_dual(self.prime_powers[1][0], self.prime_powers[1][1])
        elif mode == 4:
            poly_m[2] = sample_dual(self.prime_powers[2][0], self.prime_powers[2][1])


        return TensorMPoly(self, poly_m, mode)     

    def create_element(self, poly_m:list, mode:int):
        # assert if the element being created as the required format
        return TensorMPoly(self, poly_m, mode)

class TensorMPoly:  
    def __init__(self, field, poly_m:list, mode:int):
        assert(0 < mode < 5)
        assert(len(poly_m) == 3)

        self.mode = mode
        self.poly_m = poly_m
        self.field = field
        
         
    def __add__(self, other):
        assert(self.mode == other.mode)
        
        # bom, provavelmente apenas somas onde polinomios ou multiplos inteiros de polinomios devem ocorrer
        # nos dois elementos sendo somados  

        # como o framework é feito para um procedimento específico é bem próvavel que seja possível implementar
        # para alguma situação específica
    
        # como não temos como fazer elemento -> produto tensorial não temos como fazer o caso genérico
        return None 

    def __mul__(self, other):
        assert(abs(self.mode - other.mode) == 2)
        new_poly = [self.poly_m[0] * other.poly_m[0] % self.field.cyclo_poly[0]]
        if(self.mode == 1 and other.mode == 3):
            new_poly.append(self.poly_m[1] * other.poly_m[1] % self.field.cyclo_poly[1])
            new_poly.append(other.poly_m[2])
            new_poly = trace_mode(new_poly, self.field.prime_powers, 1)

            return self.field.create_element(new_poly, 2)

        elif(self.mode == 3 and other.mode == 1):
            new_poly.append(self.poly_m[1] * other.poly_m[1] % self.field.cyclo_poly[1])
            new_poly.append(self.poly_m[2])
            new_poly = trace_mode(new_poly, self.field.prime_powers, 1)

            return self.field.create_element(new_poly, 2)

        elif(self.mode == 2 and other.mode == 4):
            new_poly.append(other.poly_m[1])
            new_poly.append(self.poly_m[2] * other.poly_m[2] % self.field.cyclo_poly[2])
            new_poly = trace_mode(new_poly, self.field.prime_powers, 2)
            
            return self.field.create_element(new_poly, 1)

        elif(self.mode == 4 and other.mode == 2):
            new_poly.append(self.poly_m[1])
            new_poly.append(self.poly_m[2] * other.poly_m[2] % self.field.cyclo_poly[2])
            new_poly = trace_mode(new_poly, self.field.prime_powers, 2)

            return self.field.create_element(new_poly, 1)

    
    def trace(self):
        return
        
## Test ## 

Field = TesorMField([(2,2),(3,2),(5,2)])
print(Field.prime_powers)
for pair in [(2,4),(1,3)]:
    el1 = Field.random_element(pair[0])
    el2 = Field.random_element(pair[1])
    el3 = el1 * el2
    print(el3.poly_m)