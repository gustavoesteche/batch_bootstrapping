load("trace/rlwe_homo_trace.sage")
from time import perf_counter 

# Testing the slower alternative 
print("## Testing RLWE homomorphic trace ##") 

# Testing parameters
p1, n1 = 2, 1
p2, n2 = 3, 1
p3, n3 = 5, 1
n = (p1 ** n1) * (p2 ** n2) * (p3 ** n3)

B, q = 2, 10000
l = ceil(log(q, B))
q = B**l
sigma = 1
base = 2

print("Parameters")
print(p1, n1)
print(p2, n2)
print(p3, n3)
print("B, q", B, q)
print("sigma ", sigma)
print("key-switch base ", base)

gsw = GSW(n, q, sigma, B)
rlwe = RLWE(n, q, B, s=gsw.sk,sigma=sigma)

prime ,exponent = p2 , n2
#assert gsw.sk == rlwe.s

auto = find_aut(n, prime, exponent)
evk_ktensor = gsw.gen_evk_ktensor(prime, exponent, base)

f = Zx(cyclotomic_polynomial(n))
Zqx = (ZZ.quotient(q))['x']
Rb = (ZZ.quotient(B))['x'].quotient(f)

start = perf_counter()
tests = 1
for i in range(tests):
    msg = random_msg_gen(B, n)
    ciphertext = rlwe.enc(msg)
    #msg_dec = rlwe.dec(ciphertext)
    #print("mensagem: ",msg)
    #print("criptograma: ",ciphertext)
    #print("mensagem decifrada: ",msg_dec)
    #assert msg == msg_dec

    # Lets assume that the P is 1, since there are no duals
    P = 1
    evk = gsw.enc(gsw.sk)
    
    trace = homo_trace_rlwe(ciphertext, P, evk, evk_ktensor, B, q, l, n, prime, exponent, base, auto)

    trace_m = generic_trace(msg, n, prime, exponent)
    # print("traço cifrado: ",trace)
    #print("traço decifrado: ",rlwe.dec(trace))
    #print("traço da mensagem esperado: ", Zx(Rb(trace_m).lift()))
    #print("traço da mensagem puro: ", trace_m)
    assert Zx(Rb(trace_m).lift()) == rlwe.dec(trace)

end = perf_counter()
print("homomorphic trace is working as intended")
print("approximate time per comparision:  ", (end - start)/tests)

