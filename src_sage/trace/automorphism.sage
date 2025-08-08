# Testing if the automorfism 'selects' the proper elements from the field

Z = PolynomialRing(ZZ, 'x') # R1

def automorphism(poly, m, i):
    return Z(poly(x=x^i)) % Z(cyclotomic_polynomial(m))

p, n = 3, 2
poly = Z.random_element(degree = euler_phi(p**n)-1)

print("automorphisms in the galois group\n")
for i in Zmod(p**n).list_of_elements_of_multiplicative_group():
    print(i, automorphism(poly, p**n, i))

# Notice that the automorphisms that are present on the Galois group work as 
# a shuffle on the polynomial, and those there arent deform the algebraic structure

print("\nautomorphisms that arent in the galois group\n")
for i in range(0,p**n,p):
    print(i, automorphism(poly, p**n, i))