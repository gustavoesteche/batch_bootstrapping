load("./trace/trace_lib.sage")
load("./GSW/GSW.sage")
load("./CRLWE/RLWE.sage")
load("./key_switch/key_switch_lib.sage")
load("./utils/utils.sage")
load("./utils/gadget_matrix.sage")
load("./dual/dual_functions.sage")

p1, n1 = 2, 1
p2, n2 = 3, 1
p3, n3 = 5, 1
n = (p1 ** n1) * (p2 ** n2) * (p3 ** n3)

B, q = 2, 2048
l = ceil(log(q, B))
sigma = 0.1
base = 2

gsw = GSW(n, q, sigma, B)
rlwe = RLWE(n, q, B, s=gsw.sk,sigma=sigma)

prime ,exponent = p1 , n1
#assert gsw.sk == rlwe.s

for i in range(100):
    m = random_msg_gen(B, n)
    c = rlwe.enc(m)
    evk_ktensor = gsw.gen_evk_ktensor(prime, exponent, base)

    mi = prime ** exponent
    fator = n / mi
    k = 0
    for j in Zmod(mi).list_of_elements_of_multiplicative_group():
        i = j - 1
        ks_res = key_switch([rlwe.Rq(Zx(c[0].lift())(x^(i*fator+1))),rlwe.Rq(Zx(c[1].lift())(x^(i*fator+1)))], evk_ktensor[k], base, B, q, n)
        m_aut = Zx(rlwe.Rq(m(x**(i*fator+1))).lift())
        assert rlwe.dec(ks_res) == Zx(rlwe.Rq(m(x**(i*fator+1))).lift())
        k += 1

for i in range(100):
    m = random_msg_gen(B, n)
    c = gsw.dist_rlwe.sample() 
    c[1] = c[1] + m

    evk_ktensor = gsw.gen_evk_ktensor(prime, exponent, base)

    mi = prime ** exponent
    fator = n / mi
    k = 0
    for j in Zmod(mi).list_of_elements_of_multiplicative_group():
        i = j - 1
        ks_res = key_switch([rlwe.Rq(Zx(c[0].lift())(x^(i*fator+1))),rlwe.Rq(Zx(c[1].lift())(x^(i*fator+1)))], evk_ktensor[k], base, B, q, n)
        m_aut = Zx(rlwe.Rq(m(x**(i*fator+1))).lift())
        assert (ks_res[1] - ks_res[0]*rlwe.s).lift() == m_aut
        k += 1

print("OK")