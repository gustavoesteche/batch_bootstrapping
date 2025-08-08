load("./utils/gadget_matrix.sage")
load("./utils/utils.sage")
load("./CRLWE/RLWE.sage")
load("./key_switch/key_switch_lib.sage")

Zx = PolynomialRing(ZZ, "x")

'''
# data points with fixed N
print("I,N,initial_noise,noise")
tests = 200
N = (2 ** 5) * 3 * 5
for i in range(12, 25):
    noise = 0
    initial_noise = 0
    B, q = 2, 2**i
    new_rlwe = RLWE(N, q, B, sigma = 2)
    old_rlwe = RLWE(N, q, B, sigma = 2)
    for _ in range(tests):
        m = random_msg_gen(B, N)
        cypher = old_rlwe.enc(m)

        base = 2
        K = new_rlwe.enc_sk(old_rlwe.s, base)
        new_cypher = key_switch(cypher, K, base,B, q, N, new_rlwe.Rq)
        initial_noise += old_rlwe.get_noise(cypher)
        noise += new_rlwe.get_noise(new_cypher) 
        d = new_rlwe.dec(new_cypher)
        assert m == d, "The key-switch did not work"
    
    print(i,end=',')
    print(N,end=',')
    print(initial_noise/tests,end=',')
    print(noise/tests)

'''

# data points with fixed Q
print("q,N,initial_noise,noise")
tests = 100
for i in range(2, 6):
    noise = 0
    initial_noise = 0
    B, q = 2, 2**30
    N = (2 ** i) * 3 * 5
    new_rlwe = RLWE(N, q, B, sigma = 2)
    old_rlwe = RLWE(N, q, B, sigma = 2)
    for _ in range(tests):
        m = random_msg_gen(B, N)
        cypher = old_rlwe.enc(m)
        base = 2
        K = new_rlwe.enc_sk(old_rlwe.s, base)
        new_cypher = key_switch(cypher, K, base,B, q, N, new_rlwe.Rq)
        initial_noise += old_rlwe.get_noise(cypher)
        noise += new_rlwe.get_noise(new_cypher) 
        d = new_rlwe.dec(new_cypher)
        assert m == d, "The key-switch did not work"
    
    print(q,end=',')
    print(N,end=',')
    print(initial_noise/tests,end=',')
    print(noise/tests)
