load("./GSW/GSW.sage")
load("./CRLWE/RLWE.sage")
load("./trace/trace_lib.sage")
load("./trace/rlwe_homo_trace.sage")
load("./ciphertext_operations_article/ciphertext_operations.sage")

Zx = ZZ['x']

# Define parameters 
p1, n1 = 7, 1
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

rmode = 1
gmode = 3
framework = Operations(n, q, B, prime_powers, base, sigma)

Rb = (ZZ.quotient(B))['x'].quotient(cyclotomic_polynomial(n))
Rq = (ZZ.quotient(q))['x'].quotient(cyclotomic_polynomial(n))

tests = 5
depth = 10
for lk in range(tests):
    msg_rlwe = sample_msg_pkgs(B, n, prime_powers, r)
    p_msg_rlwe = pack_msg(msg_rlwe, prime_powers, n, r, rmode, Rq)

    # Packing the messages for testing
    packed_rlwe = framework.pack_rlwe(msg_rlwe, rmode)
    original_gmode = gmode
    for i in range(depth):        
        msg_gsw = sample_msg_pkgs(B, n, prime_powers, r)    
        p_msg_gsw = pack_msg(msg_gsw, prime_powers, n, r, gmode, Rq)        
        packed_gsw = framework.pack_gsw(msg_gsw, gmode) 

        if gmode % 2:
            prime, exponent = p2, n2
        else:
            prime, exponent = p3, n3

        # Testing by the trace in the packed message 
        pinverse = (prime**exponent).inverse_mod(q)
        tr = Rq(generic_trace(p_msg_rlwe * p_msg_gsw, n, prime, exponent)) * pinverse  
        tr_msg = Zx(Rb(tr).lift())

        mult = framework.ext_prod_trace(packed_rlwe, packed_gsw)
        result = framework.rlwe.dec(mult[0])
        
        pack_depth = framework.unpack_rlwe(mult)
        assert result == tr_msg, f"\n {result} \n {tr_msg}"        
        for j in range(r):
            assert Rb(msg_rlwe[j] * msg_gsw[j]) == pack_depth[j]
            msg_rlwe[j] = Zx(Rb(msg_rlwe[j] * msg_gsw[j]).lift())

        noise = framework.rlwe.get_noise(mult[0], Zx(tr.lift()))
        #print("trace on the message ", tr_msg)
        #
        
        packed_rlwe = mult
        p_msg_rlwe = result
        if gmode == 4:
            gmode = 3
        else:
            gmode = 4   
        #print(noise)
    gmode = original_gmode 
print("External product with trace is working as intended")