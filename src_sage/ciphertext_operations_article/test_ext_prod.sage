load("./GSW/GSW.sage")
load("CRLWE/RLWE.sage")
load("./trace/trace_lib.sage")
load("./trace/rlwe_homo_trace.sage")
load("./ciphertext_operations_article/ciphertext_operations.sage")

Zx = ZZ['x']

# Define parameters 
p1, n1 = 2, 5
p2, n2 = 3, 1
p3, n3 = 5, 1
prime_powers = [(p1, n1), (p2, n2), (p3, n3)]
m1 , m2, m3 = p1**n1, p2**n2, p3**n3
n = m1 * m2 * m3

r = min(euler_phi(m2), euler_phi(m3))

B, q = 2, 2**100
l = ceil(log(q, B))
q = B**l
sigma = 3.2
base = 2

framework = Operations(n, q, B, prime_powers, base, sigma)

rmode = 1
gmode = 3

f = cyclotomic_polynomial(n)
Rb = (ZZ.quotient(B))['x'].quotient(f)
Rq = (ZZ.quotient(q))['x'].quotient(f)

noise = 0
initial_noise = 0
tests = 1
depth = 1 
for _ in range(tests):
    # Sampling the messages
    msgs_rlwe = sample_msg_pkgs(B, n, prime_powers, r) # random_msg_gen(B, n)
    msgs_gsw = sample_msg_pkgs(B, n, prime_powers, r) # random_msg_gen(B, n)

    # Packing messages into the big message
    p_msgR = pack_msg(msgs_rlwe, prime_powers, n, r, rmode, Rq)
    p_msgG = pack_msg(msgs_gsw, prime_powers, n, r, gmode, Rq)

    # Packing the messages into a cipher
    packed_rlwe = framework.pack_rlwe(msgs_rlwe, rmode)
    packed_gsw = framework.pack_gsw(msgs_gsw, gmode)

    # print(gsw.get_noise(packed_gsw)) # ruido inicial gsw
    # print(rlwe.get_noise(packed_rlwe)) # ruido inicial rlwe
    c = framework.test_ext_prod(packed_rlwe, packed_gsw)
    # print("Equality: ", Zx(Rb(p_msgR * p_msgG).lift()) == rlwe.dec(c))
    # print(Zx(Rb(p_msgR * p_msgG).lift()))
    # print(rlwe.dec(c))
    
    # compute the noises
    initial_noise += framework.rlwe.get_noise(packed_rlwe[0])
    noise += (framework.rlwe.get_noise(c))

    assert Rb(p_msgR * p_msgG) == framework.rlwe.dec(c)

print("Ext-Prod working as intended")
#print("Average initial noise rlwe ", (initial_noise/tests).n())
#print("Average noise ", (noise/tests).n())