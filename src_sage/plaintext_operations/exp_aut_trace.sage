# Test experimentally if automorfisms and traces associate, it does validating the theory 

Z = PolynomialRing(ZZ, 'x')

def automorphism(f, a, cyc):
    return Z(f(x = x^a)) % cyc

def coprimes(p, m):
    coprimes_list = []
    for j in range(0, p**(m-1)):
        for i in range(1, p):
            coprimes_list.append(j * p + i)

    return coprimes_list


## teste sigma_1(sigma_2(a)) = sigma_2(sigma_1(a)) para sigma_1 e sigma_2 automorfismos :) 
'''
p, m = 3,3 ## random_prime(20), ZZ.random_element(1, 5)
cyc = Z(cyclotomic_polynomial(p ** m))
f = Z.random_element(degree = euler_phi(p ** m)-1)
aut = coprimes(p, m)
flag = True

for i in range(len(aut)):
    for j in range(i+1, len(aut)):
        if(i != j): 
            A = automorphism(automorphism(f, aut[i], cyc), aut[j], cyc)
            B = automorphism(automorphism(f, aut[j], cyc), aut[i], cyc)
            if A != B:
                print(i, j, A)
                print(j, i, B)
                print("Falso")
                flag = False

if(flag):
    print("Correto")
'''

