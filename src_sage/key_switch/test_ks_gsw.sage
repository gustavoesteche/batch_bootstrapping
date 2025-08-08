load("key_switch/key_switch_lib.sage")
load("GSW/GSW.sage")
    
Zx = ZZ['x']

def random_msg_gen(B, n):
    Rb = (ZZ.quotient(B) )['x'].quotient(cyclotomic_polynomial(n))
    return Zx((Rb.random_element()).lift())

# Setting parameters
n = 2 ** 4
B, q = 2, 2**10
l = ceil(log(q, B))
sigma = 10

# Instaciating ciphers
old_gsw = GSW(n, q, sigma, B)
new_gsw = GSW(n, q, sigma, B)

noise = 0
initial_noise = 0 

# Testing the corretude of the key-switch
tests = 200
for i in range(tests):
    m = random_msg_gen(B, n)
    c = old_gsw.enc(m)

    base = 2
    K = new_gsw.enc_sk(old_gsw.sk, base) 
    new_c = ks_gsw(c, new_gsw.Ks, K, base, new_gsw.B, new_gsw.q, new_gsw.n, new_gsw.Rq)

    initial_noise += old_gsw.get_noise(c)
    noise += new_gsw.get_noise(new_c, m)

    assert new_gsw.dec(new_c) == m

print("GSW key-switch is working as intended")
#print("Average initial noise gsw ", initial_noise/tests)
#print("Average noise after the key-switch ", noise/tests)
