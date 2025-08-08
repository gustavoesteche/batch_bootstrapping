load("./GSW/GSW.sage")
load("./CRLWE/RLWE.sage")
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

# Instanciating RGSW
gsw = GSW(n, q, sigma, B)
gmode = 3

# Instanciating RLWE
rlwe = RLWE(n, q, B, s=gsw.sk ,sigma=sigma)
rmode = 1

Rb = (ZZ.quotient(B))['x'].quotient(gsw.f) # defining the message ring for testing purposes
Rq = (ZZ.quotient(q))['x'].quotient(gsw.f)

# Sampling messages
msgs_rlwe = sample_msg_pkgs(B, n, prime_powers, r)
msgs_gsw = sample_msg_pkgs(B, n, prime_powers, r)

# Packing the messages for testing
p_msgR = pack_msg(msgs_rlwe, prime_powers, n, r, rmode, Rq)
p_msgG = pack_msg(msgs_gsw, prime_powers, n, r, gmode, Rq)

res = Rb(0)
for i in range(r):
    res = res + Rb(msgs_rlwe[i]) * Rb(msgs_gsw[i]) * Rb(Zx(x**(i * (p1**n1) * (p2**n2))))

# Testing by the trace in the packed message 
tr = generic_trace(p_msgR * p_msgG, n, p2, n2) 
tr_msg = Zx(Rb(tr).lift())

assert tr_msg == res
print("Trace over packed messages is working as intended") 