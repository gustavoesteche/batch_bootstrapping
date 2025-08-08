load("./GSW/GSW.sage")
load("./dual/dbasis.sage")
load("./ciphertext_operations_article/ciphertext_operations.sage")

# Testing 
p1, n1 = 2, 5
p2, n2 = 3, 1
p3, n3 = 5, 1

prime_powers = [(p1, n1), (p2, n2), (p3, n3)]
m1 , m2, m3 = p1**n1, p2**n2, p3**n3
n = m1 * m2 * m3

r = min(euler_phi(m2), euler_phi(m3))

B, q = 2, 2**100
l = ceil(log(q, B))
sigma = 3.2
base = 2

framework = Operations(n, q, B, prime_powers, base, sigma)
mode = 3

Rb = ( ZZ.quotient(B) )['x'].quotient(cyclotomic_polynomial(n))
Rq = ( ZZ.quotient(q) )['x'].quotient(cyclotomic_polynomial(n))

packages = sample_msg_pkgs(B, n, prime_powers, r)
#print("packages ", packages,'\n')

packed = framework.pack_gsw(packages, mode)
#print("pack ", packed,'\n') 

#print("error of the calculated message ", gsw.get_noise(packed))

#unpacked = framework.unpack_gsw(packed)
#print("unpacked ", unpacked)

for i in range(r):
    msg_exp = Rb(packages[i])
    #print(msg_exp)
    #print(unpacked[i])
    #assert msg_exp == Rb(unpacked[i]), "{} \n {}".format(msg_exp, unpacked[i])

print("GSW pack-unpack working as intended")