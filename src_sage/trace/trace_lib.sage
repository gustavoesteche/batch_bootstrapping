Zx = PolynomialRing(ZZ, 'x')

# Compute the trace of p^m to p^(m-1)
def trace_pm_to_pm_1(p:int, m:int, poly):
    '''Computates the trace from the integer ring Z[x]/phi_p^m(x) to
    the integer ring Z[x]/phi_p^(m-1)(x)'''    
    f = Zx(cyclotomic_polynomial(p^m))
    response = Zx(0)
    
    # If m == 1,this only defines p-1 automorphisms
    if(m == 1):
        # Sum all the exponents that are coprime with p, [1,p-1] 
        for i in range(1,p):
            response = response + Zx(poly(x=x^(i))) % f 
    
    else:    
        # Sum using the defined automorphism 
        for i in range(p):
            response = (response + Zx(poly(x=x^(i*(p^(m-1)) + 1))) % f) % f 

        # Finishing the automorphism to contain only elements in the proposed field 
        response = response(x=x^(1/p))

    return Zx(response)

# Compute the trace of p^m to p^1, using the tower extension property
# it computates the trace using trace_pm_to_pm_1 
def trace_pm_to_1(p, m, poly):
    '''Computates the trace from the integer ring Z[x]/phi_p^m(x) to
    the integers'''   
    # print("original polynomial: ",poly, "\n")
    for i in range(m, 0, -1):
        poly = trace_pm_to_pm_1(p, i, poly)
        # print("trace from {}^{} to {}^{}:".format(p,i,p, i-1),poly,"\n")
    
    return poly

# A fast way to implement the trace, instead of computing the sum of all
# automorphisms, it already uses the theoretical result
def fast_trace_pm_to_pm_1(p, m, poly):
    '''Computates the trace more efficiently from the integer ring 
    Z[x]/phi_p^m(x) to the integer ring Z[x]/phi_p^(m-1)(x)''' 
    f = Zx(cyclotomic_polynomial(p ** m) )

    response = Zx(0)
    if(m == 1):
        for i in range(1,p):
            response = response + Zx(poly(x=x^(i))) % f 
            
    else:
        for i in range(0, poly.degree()+1, p):
            response = Zx(response + Zx(poly[i] * p * x^i)) % f

        response = response(x=x^(1/p))
    return Zx(response) 

# Compute the trace of p^m to p^1, using the tower extension property
# it computates the trace using fast_trRace_pm_to_pm_1 
def fast_trace_pm_to_1(p, m, poly):
    '''Computates the trace more efficiently from the integer
    ring Z[x]/phi_p^m(x) to the integers'''   
    response = poly
    for i in range(m, 0, -1):
        response = fast_trace_pm_to_pm_1(p, i, response)

    return Zx(response)

# The trivial form to compute the trace for to the integers
def trivial_trace_m_to_1(poly, m):
    '''Computates the trace from the integer ring Z[x]/phi_m(x) to the integers, 
    summing all the automorphisms defined by the exponents which are coprime to m  
    '''
    f = cyclotomic_polynomial(m)
    coprimes = [x for x in range(1, m) if gcd(x, m) == 1]

    response = 0
    for i in coprimes:
        term = Zx(poly(x=x^i)) % f
        response = response + term
    return Zx(response)

def generic_trace_vchange(poly, m, p, n):
    # function that computes trace from K_m -> K_m/(p**n) with the variable change 
    m_i = p ** n 
    factor = m / m_i
    res = Zx(0)
    for i in Zmod(m_i).list_of_elements_of_multiplicative_group():
       res = res + aut(poly, m, i*factor + 1)
       
    res = Zx(res(x=x^(m_i.inverse_mod(factor)))) % Zx(cyclotomic_polynomial(factor))
    return res      

def find_aut(m, p, n):
    mi = p ** n
    fator = m/mi
    
    k = p-(fator % p).inverse_mod(p)
    rest = [i for i in range(p) if i!=k]
    p_power = mi//p
    aut = [(i*p + j)*fator+1 for i in range(p_power) for j in rest]
    return aut

# Very inefficient approach it will work for now though if I need to test something
def create_table(m, p, n):
    mi = p ** n
    fator = m//mi
    table = vector(Zx, [0] * euler_phi(fator))
    for i in range(euler_phi(fator)):
        table[i] = Zx(x**(i*mi)) % Zx(cyclotomic_polynomial(m))

    return table
    
def reverse_aut(poly, table):
    res = Zx(0)
    for i in range(len(table)-1,-1,-1):
        cf = Zx(poly // table[i])[0]
        res = res + cf * Zx(x**i)
        poly = poly - cf * table[i]

    return res

def generic_trace(poly, m, p, n):
    f = Zx(cyclotomic_polynomial(m))
    auto = find_aut(m, p, n) 
    ans = Zx(0)
    for i in auto:
        ans = ans + (big_aut(poly, i, Zx) % f)
    
    return ans   

def generic_trace_tower(poly, m, p, n):
    m_i = p ** n 
    f = Zx(cyclotomic_polynomial(m))
    while(m % p == 0):
        ans = Zx(0)
        auto = find_aut_tower(m, p)
        for i in auto:
            ans = ans + (big_aut(poly, i, Zx) % f)
        m = m/p
        n = n - 1
        poly = ans 
    
    return ans    

'''
# Testing Example

from time import time

p, m = 5, 3
polynomial = R.random_element(degree = euler_phi(p ** m) - 1)

print("sample polynomial: ", polynomial)

start = time()
print("trivial trace computation", trivial_trace_m_to_1(polynomial, p ** m))
end = time()
print("trivial time ", end - start) 

start = time()
print("summing all automorphisms in tower method", trace_pm_to_1(p, m, polynomial))
end = time()
print("normal tower time ", end-start) 

start = time()
print("quick trace tower method", fast_trace_pm_to_1(p,m,polynomial))
end = time()
print("power-up tower time ", end-start) 
'''