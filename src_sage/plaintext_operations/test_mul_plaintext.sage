
load("./plaintext_operations/plaintext_operations.sage")
load("./trace/trace_lib.sage")
load("./dual/dbasis.sage")

Zx = PolynomialRing(ZZ, 'x') # R1
Qx = PolynomialRing(QQ, 'x') # K1

R.<x, y, z> = PolynomialRing(ZZ) # R
K.<x, y, z> = PolynomialRing(QQ) # K

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

def compose_poly(prime_powers:list, decomposed_poly:list, m:int):
    poly = 1
    for i in range(len(decomposed_poly)):
        poly = poly * decomposed_poly[i](x = x^(m/(prime_powers[i][0] ** prime_powers[i][1])))

    return Zx(poly) % Zx(cyclotomic_polynomial(m))

# Test 
p1, n1 = 5,1
p2, n2 = 3,2
p3, n3 = 2,2

m1 = p1 ** n1
m2 = p2 ** n2
m3 = p3 ** n3

m = m1 * m2 * m3
r = min(euler_phi(m2), euler_phi(m3)) # amount of elements packed
prime_powers = [(p1,n1),(p2,n2),(p3,n3)]

f1 = R(cyclotomic_polynomial(m1))(x=x)
f2 = R(cyclotomic_polynomial(m2))(x=y)
f3 = R(cyclotomic_polynomial(m3))(x=z)

db_x = canon_dbasis(p2, n2)
db_y = []
for i in range(len(db_x)):
    db_y.append(K(db_x[i](x=y)))

K_ = FieldK(prime_powers)

epochs = 10
for _ in range(epochs):
    # print("#",end="")  
    pol1 = K_.random_element(1)  # without loss of generality it is possible to fix the modes of the plaintext 
    pol2 = K_.random_element(3)  # since the multiplcation between 1/3 and 2/4 are congruent.
    pol3 = pol1 * pol2

    #print("element 1: ", pol1.polys, "mode", pol1.mode)
    #print("element 2: ",pol2.polys, "mode", pol2.mode)
    #print("pol1 * pol2: ",pol3.polys)

    answer = K_.compose_poly(pol3)
    #print("\ncomposto em m1*m3: ",answer,"\n")

    p1_cp = K(0)
    p2_cp = K(0)

    # unpacking the elements 
    for i in range(r):
        p1_cp = (p1_cp + K(pol1.polys[i]) * K(y^i))

    for i in range(r):
        p2_cp = p2_cp + K(pol2.polys[i]) * db_y[i] * K(z^i)  

    # multiplying both elements
    p2_cp = R(m2 * p2_cp) # mapping out the elements from the dual
    p_cp = R(p1_cp * p2_cp) # multiplying composed polynomials 
    p_cp = ((p_cp.mod(f1)).mod(f2)).mod(f3)
    #print("tensor outer product: ",p_cp)

    # composing the element
    p_cp = Zx(p_cp(x=x^(m2*m3), y=x^(m1*m3), z=x^(m1*m2))) % Zx(cyclotomic_polynomial(m)) 
    #print("after the isomorphism: ", p_cp)

    # trace computation and division
    trace = generic_trace(p_cp, m, p2,n2)/m2
    # print("final trace: ", trace)
    assert trace == answer, str(trace)+""+str(answer)

print("Multiplication of plaintexts is working as intended :)\n")
