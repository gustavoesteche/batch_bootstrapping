# The real use for this implementation is to understand the varialble change to be
# implemented

# TODO: needs to correction

load("./trace/trace_lib.sage")

def generic_tower_trace(poly, m, p):
    factor = m / p
    res = Zx(0)
    for i in range(p):
       res = res + aut(poly, m, i*factor + 1)

    if (m / p) % p == 0:
        res = Zx(res(x=x**(1/p))) % Zx(cyclotomic_polynomial(factor))
    else:
        res = Zx(res(x=x^(p.inverse_mod(factor)))) % Zx(cyclotomic_polynomial(factor))

    return res      

def hello(poly, m, p, n):
    for i in range(n):
        poly = generic_tower_trace(poly, m, p)
        m = m / p
    return poly

p1, n1 = 2, 1
p2, n2 = 3, 1
p3, n3 = 5, 1

m1 = p1 ** n1
m2 = p2 ** n2
m3 = p3 ** n3

m = m1 * m2 * m3

poly = Zx.random_element(degree=euler_phi(m)-1)

print(generic_trace(poly, m, p1, n1))
print(hello(poly, m, p1, n1))