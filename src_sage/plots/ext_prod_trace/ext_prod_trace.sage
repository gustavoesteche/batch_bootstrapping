load("./GSW/GSW.sage")
load("./CRLWE/RLWE.sage")
load("./trace/trace_lib.sage")
load("./trace/rlwe_homo_trace.sage")
load("./ciphertext_operations_article/ciphertext_operations.sage")

Zx = ZZ['x']

# Setting parameters 
p1, n1 = 2, 3
p2, n2 = 13, 1
p3, n3 = 3, 2  
prime_powers = [(p1, n1), (p2, n2), (p3, n3)]
m1 , m2, m3 = p1**n1, p2**n2, p3**n3
n = m1 * m2 * m3

r = min(euler_phi(m2), euler_phi(m3))

B, q = 2, 2**50
l = ceil(log(q, B))
q = B**l
sigma = 3.2
base = 2

rmode = 1
gmode = 3
framework = Operations(n, q, B, prime_powers, base, sigma)
print("Instância criada")

Rb = (ZZ.quotient(B))['x'].quotient(cyclotomic_polynomial(n)) 
Rq = (ZZ.quotient(q))['x'].quotient(cyclotomic_polynomial(n)) 

# print("noise,depth")

tests = 1
depth = 5
#noise = [0] * depth
for _ in range(tests):
    msg_rlwe = sample_msg_pkgs_depth(B, n, prime_powers, r)
    p_msg_rlwe = pack_msg(msg_rlwe, prime_powers, n, r, rmode, Rq)  
    packed_rlwe = framework.pack_rlwe(msg_rlwe, rmode)

    msgs_gsw = [] 
    for i in range(depth):
        msg = sample_msg_pkgs_depth(B, n, prime_powers, r)
        msgs_gsw.append(msg)

    o_gmode = gmode
    for i in range(depth):        
        #print("profundidade ", i)
        if gmode % 2:
            prime, exponent = p2, n2
            fator = n / (p3**n3)
        else:
            prime, exponent = p3, n3
            fator = n / (p2**n2)

        # Testing by the trace in the packed message
        p_msg_gsw = pack_msg(msgs_gsw[i], prime_powers, n, r, gmode, Rq)
        packed_gsw = framework.pack_gsw(msgs_gsw[i], gmode) 
        tr = Rq(generic_trace(p_msg_rlwe * p_msg_gsw , n, prime, exponent)) * (prime**(exponent)).inverse_mod(q) 
        tr = Zx(tr.lift())
        tr_msg = Zx(Rb(tr).lift())

        mult = framework.ext_prod_trace(packed_rlwe, packed_gsw)
        result = framework.rlwe.dec(mult[0])

        #noise[i] = max(noise[i], framework.rlwe.get_noise(mult[0], tr))
        print(framework.rlwe.get_noise(mult[0], tr), end=',')
        assert result == tr_msg

        packed_rlwe = mult
        p_msg_rlwe = tr
        if gmode == 4:
            gmode = 3
        else:
            gmode = 4   
    gmode = o_gmode
    
#for i in range(depth):
    #print((noise[i]).n(),end=",")    
    #print`(i+1)
