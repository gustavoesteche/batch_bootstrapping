Zx = PolynomialRing(ZZ, 'x')
Qx = PolynomialRing(QQ, 'x') 

def automorphism(poly, m, i):
    return Zx(poly(x=x^i)) % Zx(cyclotomic_polynomial(m))

# function that computes trace from K_m -> K_m/(p**n)
def generic_trace(poly, m, p, n):
    m_i = p ** n 
    factor = m / m_i
    res = Zx(0)
    for i in Zmod(m_i).list_of_elements_of_multiplicative_group():
       res = res + automorphism(poly, m, i*factor + 1)

    res = Zx(res(x=x^(m_i.inverse_mod(factor)))) % Zx(cyclotomic_polynomial(factor))
    return res  

# The trivial form to compute the trace for to the integers
def trivial_trace_m_to_1(poly, m):
    '''Computates the trace from the integer ring Z[x]/phi_m(x) to the integers, 
    summing all the automorphisms defined by the exponents which are coprime to m  
    '''
    f = cyclotomic_polynomial(m)
    coprimes = [x for x in range(1, m) if gcd(x, m) == 1]

    response = 0
    for i in coprimes:
        term = Qx(poly(x=x^i)) % f
        response = response + term
    return Qx(response)

def separate_multiples(poly, k):
    poly_l = list(poly)
    ans = Zx(0)

    for i in range(0,len(poly_l),k):
        ans = ans + poly_l[i] * Zx(x^i)
    return ans

# Test if the tower of generic traces leads to the correct result
# While also testing each step with the previous knowlegde of generic trace

p1, n1 = 5, 1
p2, n2 = 3, 1
p3, n3 = 2, 1

m1 = p1 ** n1
m2 = p2 ** n2
m3 = p3 ** n3

m = m1 * m2 * m3

e = 100
for i in range(e):
    poly1 = Zx.random_element(degree = euler_phi(m1)-1)
    poly2 = Zx.random_element(degree = euler_phi(m2)-1)
    poly3 = Zx.random_element(degree = euler_phi(m3)-1)

    comp = Zx(poly1(x=x^(m2*m3))) * Zx(poly2(x=x^(m1*m3))) * Zx(poly3(x=x^(m1*m2))) % Zx(cyclotomic_polynomial(m))

    tota_trace = trivial_trace_m_to_1(poly1, m1)*trivial_trace_m_to_1(poly2, m2) * trivial_trace_m_to_1(poly3, m3)
    
    trace_123_12 = generic_trace(comp, m, p3, n3)
    #print("Trace R123 -> R12",trace_123_12)
    
    t123_12 = Zx(trivial_trace_m_to_1(poly3, m3)* poly1(x=x^(m2))) * Zx(poly2(x=x^(m1))) % Zx(cyclotomic_polynomial(m1*m2))
    assert trace_123_12 == t123_12
    

    trace_12_1 = generic_trace(trace_123_12, m1*m2, p2, n2)
    #print("Trace R12 -> R1",trace_12_1)
    t12_1 = poly1 * trivial_trace_m_to_1(poly2, m2) * trivial_trace_m_to_1(poly3, m3)
    assert trace_12_1 == t12_1
    
    trace = trivial_trace_m_to_1(trace_12_1, m1)

    #print("Correct trace computation: ",tota_trace)
    #print("Generic trace tower trace computation: ",trace)
    assert tota_trace == trace
    #print("Test {} passed".format(i+1))

print("generic trace working as intended")