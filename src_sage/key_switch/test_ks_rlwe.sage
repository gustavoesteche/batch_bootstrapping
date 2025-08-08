load("./key_switch/key_switch_lib.sage")
load("./CRLWE/RLWE.sage")

Zx = PolynomialRing(ZZ, "x")

def random_msg_gen(B, n):
    Rb = (ZZ.quotient(B) )['x'].quotient(cyclotomic_polynomial(n))
    return Zx((Rb.random_element()).lift())

noise = 0
initial_noise = 0

# Setting parameters
N = (2 ** 5) * 3 * 5
B, q = 2, 2**11
base = 2

Rb = (ZZ.quotient(B) )['x'].quotient(cyclotomic_polynomial(N))
Rq = (ZZ.quotient(q) )['x'].quotient(cyclotomic_polynomial(N))

# Testing the corretude of the key-switch
tests = 200
for _ in range(tests):
    m = random_msg_gen(B, N) # Sampling message

    new_rlwe = RLWE(N, q, B, sigma = 2)
    old_rlwe = RLWE(N, q, B, sigma = 2)

    c = old_rlwe.enc(m)

    base = 2
    K = new_rlwe.enc_sk(old_rlwe.s, base)
    new_c = ks_rlwe(c, K, base, B, q, N, Rq)

    initial_noise += old_rlwe.get_noise(c)
    noise += new_rlwe.get_noise(new_c) 

    d = new_rlwe.dec(new_c)
    assert m == d

print("RLWE key-switch is working as intended")
#print("Average initial noise rlwe ", initial_noise/tests)
#print("Average noise ", noise/tests)