load("./plaintext_operations/plaintext_operations.sage")
load("./trace/trace_lib.sage")
load("./dual/dbasis.sage")

Zx = PolynomialRing(ZZ, 'x') # R1
Qx = PolynomialRing(QQ, 'x') # K1

R.<x, y, z> = PolynomialRing(ZZ) # R
K.<x, y, z> = PolynomialRing(QQ) # K


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

f1 = Zx(cyclotomic_polynomial(m1))(x=x)
f2 = Zx(cyclotomic_polynomial(m2))(x=y)
f3 = Zx(cyclotomic_polynomial(m3))(x=z)

K_ = FieldK(prime_powers)

epochs = 10 # this test needs a fix
for _ in range(epochs):
    # print("#",end="")  
    pol1 = K_.random_element(1)  
    pol2 = K_.random_element(1)  
    pol3 = pol1 + pol2

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
        p2_cp = p2_cp + K(pol2.polys[i]) * K(y^i) 

    # adding both elements
    p_cp = R(p1_cp + p2_cp) # adding composed polynomials
    p_cp = ((p_cp.mod(f1)).mod(f2)).mod(f3)
    #print("tensor outer product: ",p_cp)

    # composing the element to R1xR2
    p_cp = Z(p_cp(x=x^(m2), y=x^(m1))) % Z(cyclotomic_polynomial(m1*m2)) 
    #print("after the isomorphism: ", p_cp)

    #print(answer)
    #print(p_cp)
    assert p_cp == answer 

print("Addition of plaintexts is working as intended :)\n")