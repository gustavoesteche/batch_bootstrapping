load("trace/rlwe_homo_trace.sage")
from time import perf_counter 

# Testing the faster alternative
print("## Testing RLWE homomorphic trace tower ##") 

# Test parameters
p1, n1 = 1, 1
p2, n2 = 3, 1
p3, n3 = 5, 1
n = (p1 ** n1) * (p2 ** n2) * (p3 ** n3)

B, q = 2, 2**14
l = ceil(log(q, B))
q = B**l
sigma = 1
base = 2

print("Parameters")
print(p1, n1)
print(p2, n2)
print(p3, n3)

print("B =",B, "q =", q)
print("sigma ", sigma)
print("key-switch base ", base)

gsw = GSW(n, q, sigma, B)
rlwe = RLWE(n, q, B, s=gsw.sk,sigma=sigma)

prime ,exponent = p2 , n2
#assert gsw.sk == rlwe.s

auto = construct_aut_tower(n, prime)
evk_ktensor = gsw.gen_evk_ktensor_tower(prime, exponent, base)

noise = []
tests = 1
start = perf_counter()
for i in range(tests):
    P = 1  # not using duals for testing the function
    msg = Zx([-1,1])#random_msg_gen(B, n)
    print(msg)
    ciphertext = vector(rlwe.enc(msg))*P    
    evk = gsw.enc(gsw.sk)

    trace = fast_homo_trace_rlwe(ciphertext, P, evk, evk_ktensor, B, q, l, n, prime, exponent, base, auto)
    noise.append(rlwe.get_noise(trace))

    f = Zx(cyclotomic_polynomial(n))
    Rb = (ZZ.quotient(B))['x'].quotient(f)
    
    trace_m = generic_trace(msg*P, n, prime, exponent) 
    #print("traço da mensagem esperado na decriptação : ", Zx(Rb(trace_m).lift()))
    #print("traço da mensagem puro: ", trace_m)
    #print(rlwe.dec(trace))
    #print(Zx(Rb(trace_m).lift()))
    assert Zx(Rb(trace_m).lift()) == rlwe.dec(trace) 

end = perf_counter()
print("homomorphic tower trace is working as intended")
print("approximate time per comparision:  ", (end - start)/tests)
print("Average Normalized Noise at the end of the trace ", sum(noise)/tests)