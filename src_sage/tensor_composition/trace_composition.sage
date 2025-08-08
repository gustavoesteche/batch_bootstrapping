load("./trace/trace_lib.sage")

Zx = PolynomialRing(ZZ, 'x')

# Compose the polynomial in the original polynomial ring,
# performing the isomorfism and then the multiplication(ye
def compose_poly(factors_m:list, poly_m:list, m:int):
    poly = 1
    for i in range(len(poly_m)):
        poly = poly * poly_m[i](x^(m/(factors_m[i][0] ** factors_m[i][1])))

    return Zx(poly) % Zx(cyclotomic_polynomial(m))

# Compute the trace using the decomposed polynomials 
def composed_trace(factors_m:list, poly_m:list):
    trace = 1
    for i in range(len(factors_m)):
        trace = trace * trivial_trace_m_to_1(poly_m[i], factors_m[i][0] ** factors_m[i][1])
    return trace

'''
# Testing example

poly_m = []
m1 = 3
m2 = 5
m = m1* m2

poly_m.append(Zx.random_element(degree = euler_phi(m1)))
poly_m.append(Zx.random_element(degree = euler_phi(m2))) 

factors_m = list(factor(m))

print("Polinomios nos aneis de potencias de primos", poly_m)
print("Fatorizacao de m em potencias de primos",factors_m)
print("Traço composto", composed_trace(factors_m, poly_m))
composed_poly = compose_poly(factors_m, poly_m, m)
print("traço trivial", trivial_trace_m_to_1(composed_poly, m))
'''