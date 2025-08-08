load("trace/trace_lib.sage")
Zx = ZZ['x']

def find_aut(m, p, n):
    mi = p ** n
    fator = m/mi
    
    k = p-(fator % p).inverse_mod(p)
    rest = [i for i in range(p) if i!=k]
    aut = [(i*p + j)*fator+1 for i in range(p**(n-1)) for j in rest]
    return aut

'''
# testing automorphisms
prime_powers = [(3,1),(2,2),(5,1)]

n = 1
for i in prime_powers:
    n *= i[0] ** i[1]

f = Zx(cyclotomic_polynomial(n))
print(f)

poly = Zx.random_element(degree = euler_phi(n)-1) % f

print(poly)
# test 1
for index in prime_powers:
    print("#############", index,"#############")
    auto = find_aut(n, index[0], index[1])
    for i in auto:
        #ok = Zx(f(x**(i)))
        ok = Zx(poly(x**(i)))
        print(ok % f)
'''

# testing old trace function
def generic_trace_old(poly, m, p, n):
    m_i = p ** n 
    factor = m / m_i
    res = Zx(0)
    for i in Zmod(m_i).list_of_elements_of_multiplicative_group():
        res = res + aut(poly, m, i*factor + 1)

    res = Zx(res(x=x^(m_i.inverse_mod(factor)))) % Zx(cyclotomic_polynomial(factor))
    return res      

prime_powers = [(3,1),(2,1),(5,1)]
p, n = 2, 1

m = 1
for i in prime_powers:
    m *= i[0] ** i[1]

f = Zx(cyclotomic_polynomial(m))

poly = Zx.random_element(degree = euler_phi(m)-1) % f

print("poly ", poly)
print(generic_trace_old(poly, m, p, n))
print(generic_trace_new(poly, m, p, n))

