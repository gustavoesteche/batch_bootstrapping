load("./GSW/GSW.sage")
load("CRLWE/RLWE.sage")

Zx = ZZ['x']


n = 2 ** 4 * 3 ** 2
B, q = 2, 2**40
l = ceil(log(q, B))
q = B**l
sigma = 2

def random_msg_gen(B, n):
    k = ZZ.random_element(0, euler_phi(n))
    k1 = ZZ.random_element(1,2)
    if k == euler_phi(n):
        return 0
    else:
        if k1 == 1:
            return Zx(x^k)
        else:
            return Zx(-x^k)

# Instanciating RGSW
gsw = GSW(n, q, sigma, B)
gmode = 3

# Instanciating RLWE
rlwe = RLWE(n, q, B, s=gsw.sk ,sigma=sigma)
rmode = 1

Rb = (ZZ.quotient(B))['x'].quotient(gsw.f)

noise = 0
initial_noise = 0
depth = 30
m = random_msg_gen(B, n)
c = rlwe.enc(m)
e = []
for _ in range(depth):
    m1 = random_msg_gen(B, n)
    c1 = gsw.enc(m1)

    m = m * m1
    c = gsw.ext_prod(c, c1)
    print(rlwe.dec(c) == Zx(Rb(m).lift()))
 
    e.append(rlwe.get_noise(c, m % Zx(cyclotomic_polynomial(n))))

print(e)
