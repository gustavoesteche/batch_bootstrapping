load("utils/utils.sage")
load("utils/RLWE_distribution.sage")
load("utils/gadget_matrix.sage")
load("GSW/GSW.sage")
load("key_switch/key_switch_lib.sage")

Zx = ZZ['x']

# Defining parameters

# data points with fixed N
print("I,N,initial_noise,noise")
for i in range(10,25):
    n = 2 ** 4
    B, q = 2, 2**i
    l = ceil(log(q, B))
    sigma = 10
    old_gsw = GSW(n, q, sigma, B)
    new_gsw = GSW(n, q, sigma, B)
    tests = 100
    noise = 0
    initial_noise = 0 
    for _ in range(tests):
        m = random_msg_gen(B, n)
        c = old_gsw.enc(m)
        base = 2
        K = new_gsw.enc_sk(old_gsw.sk, base) 
        new_c = key_switch_gsw(new_gsw, K, base, c)
        initial_noise += old_gsw.get_noise(c)
        noise += new_gsw.get_noise(new_c, m)
        assert new_gsw.dec(new_c) == m

    print(i,end=',')
    print(N,end=',')
    print(initial_noise/tests,end='')
    print(noise/tests)

'''
# data points with fixed Q
print("q,N,initial_noise,noise")
for i in range(10,25):
    n = 2 ** i
    B, q = 2, 2**20
    l = ceil(log(q, B))
    sigma = 10
    old_gsw = GSW(n, q, sigma, B)
    new_gsw = GSW(n, q, sigma, B)
    tests = 100
    noise = 0
    initial_noise = 0 
    for _ in range(tests):
        m = random_msg_gen(B, n)
        c = old_gsw.enc(m)
        base = 2
        K = new_gsw.enc_sk(old_gsw.sk, base) 
        new_c = key_switch_gsw(new_gsw, K, base, c)
        initial_noise += old_gsw.get_noise(c)
        noise += new_gsw.get_noise(new_c, m)
        assert new_gsw.dec(new_c) == m

    print(i,end=',')
    print(N,end=',')
    print(initial_noise/tests,end='')
    print(noise/tests)'''