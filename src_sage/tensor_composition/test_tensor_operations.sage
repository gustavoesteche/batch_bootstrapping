load("./trace/trace_lib.sage")
load("./tensor_composition/trace_composition.sage")

# Testing the properties of the binary operations gained by the tensor
# decomposition and also the trace property
Zx = PolynomialRing(ZZ, 'x')

# The decompose operation is not attainable
class test_tensor_operations():
    def __init__(self, m1 = None, m2 = None):
        if(m1 == None and m2 == None):
            self.m1 = next_prime(ZZ.random_element(2, 5)) ** ZZ.random_element(1, 3)
            while(True):
                self.m2 = next_prime(ZZ.random_element(2, 5)) ** ZZ.random_element(1, 3)
                if(self.m2 % self.m1 != 0 and self.m1 % self.m2 != 0):
                    break 
        else:
            self.m1 = m1
            self.m2 = m2

        # print("test with the first prime power: ", self.m1," and the second prime power: ", self.m2)
        self.m = self.m1 * self.m2
        self.poly_m = []
        self.poly_m.append(Zx.random_element(degree = euler_phi(self.m1)-1))
        self.poly_m.append(Zx.random_element(degree = euler_phi(self.m2)-1))
        self.factors_m = list(factor(self.m))
        self.composed_poly = Zx(compose_poly(self.factors_m, self.poly_m, self.m))

    def test_trace(self):
        # when testing the case make sure to not use large values because of the trivial trace function
        trivial_trace = trivial_trace_m_to_1(self.composed_poly, self.m)
        comp_tr = composed_trace(self.factors_m, self.poly_m)
        assert (trivial_trace == comp_tr)

    # a_1 tensor b + a_2 tensor b = (a_1 + a_2) tensor b
    def test_1(self):
        other_poly = []
        other_poly.append(Zx.random_element(degree = euler_phi(self.m1)-1))
        other_poly.append(self.poly_m[1])
        
        poly_1 = Zx(compose_poly(self.factors_m,other_poly,self.m)) + Zx(compose_poly(self.factors_m,self.poly_m,self.m))
        other_poly[0] = Zx(other_poly[0]) + Zx(self.poly_m[0])
        poly_2 = Zx(compose_poly(self.factors_m,other_poly,self.m))
        assert (poly_1 == poly_2)

    # a tensor b_1 + a tensor b_2 = a tensor (b_1 + b_2)
    def test_2(self):
        other_poly = []
        other_poly.append(self.poly_m[0])
        other_poly.append(Zx.random_element(degree = euler_phi(self.m2)-1))
        
        poly_1 = Zx(compose_poly(self.factors_m,other_poly,self.m)) + Zx(compose_poly(self.factors_m,self.poly_m,self.m))
        other_poly[1] = Zx(other_poly[1]) + Zx(self.poly_m[1])
        poly_2 = Zx(compose_poly(self.factors_m,other_poly,self.m))
        assert (poly_1 == poly_2)

    # e(a tensor b) = (ea) tensor b = a tensor (eb) 
    def test_3(self):
        e = ZZ.random_element(1, 51)
        other_poly_1 = deepcopy(self.poly_m)
        other_poly_2 = deepcopy(self.poly_m)

        other_poly_1[0] = other_poly_1[0] * e
        poly_1 = Zx(compose_poly(self.factors_m,other_poly_1,self.m))

        other_poly_2[1] = other_poly_2[1] * e
        poly_2 = Zx(compose_poly(self.factors_m,other_poly_2,self.m))

        assert (poly_1 == poly_2 == e * self.composed_poly)     

    # (a_1 tenor b_1)(a_2 tenor b_2) = (a_1 a_2) tensor (b_1 b_2) 
    def test_4(self):
        p1 = []
        p1.append(Zx.random_element(degree = euler_phi(self.m1)-1))
        p1.append(Zx.random_element(degree = euler_phi(self.m2)-1))
        
        p1_c = Zx(compose_poly(self.factors_m,p1,self.m))
        poly_left = (p1_c * self.composed_poly) % (cyclotomic_polynomial(self.m))
        
        p1[0] = (p1[0] * self.poly_m[0]) 
        p1[1] = (p1[1] * self.poly_m[1])
        poly_right = Zx(compose_poly(self.factors_m,p1,self.m))
        assert (poly_left == poly_right), "left polynomial {} \n\n right polynomial {}".format(poly_left, poly_right) 
        

    def run_tests(self):
        self.test_1()
        self.test_2()
        self.test_3()
        self.test_4()
        self.test_trace()
        #print("Tests passed")

epochs = 1
for i in range(epochs):
    Test = test_tensor_operations()
    Test.run_tests() 

print("tensor operations are working as intended")