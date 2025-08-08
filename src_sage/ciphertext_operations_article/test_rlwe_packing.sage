load("./ciphertext_operations_article/ciphertext_operations.sage")

# Testing 
p1, n1 = 2, 5
p2, n2 = 3, 1
p3, n3 = 5, 1

prime_powers = [(p1, n1), (p2, n2), (p3, n3)]
m1 , m2, m3 = p1**n1, p2**n2, p3**n3
n = m1 * m2 * m3

r = min(euler_phi(m2), euler_phi(m3))

B, q = 2, 2**50
l = ceil(log(q, B))
sigma = 3.2
base = 2

Rb = ( ZZ.quotient(B) )['x'].quotient(cyclotomic_polynomial(n))
Rq = ( ZZ.quotient(B) )['x'].quotient(cyclotomic_polynomial(n))

framework = Operations(n, q, B, prime_powers, base, sigma)

for test_mode in [1,2,3,4]:
    packages = sample_msg_pkgs(B, n, prime_powers, r)
    #print("packages ", packages,'\n')

    packed = framework.pack_rlwe(packages, test_mode)
    #print("pack ", packed,'\n') 

    unpacked = framework.unpack_rlwe(packed)
    #print("unpacked ", unpacked)

    for i in range(r):
        msg_exp = Zx(Rb(packages[i]).lift())
        #print(unpacked[i])
        #print(msg_exp)
        #print(msg_exp == unpacked[i])
        assert msg_exp  == unpacked[i], "message: {} \n result; {} \n operating mode: {} ".format(msg_exp, unpacked[i], test_mode)

print("RLWE pack-unpack working as intended")
