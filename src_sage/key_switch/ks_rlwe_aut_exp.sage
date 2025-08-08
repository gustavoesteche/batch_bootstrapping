load("utils/utils.sage")
load("utils/RLWE_distribution.sage")
load("utils/gadget_matrix.sage")
load("CRLWE/RLWE.sage")
load("key_switch/RLWE_KS.sage")
load("trace/trace_lib.sage")

Zx = ZZ['x']

# Defining parameters
p1, n1 = 2, 2
p2, n2 = 3, 2
p3, n3 = 7, 2
n = (p1 ** n1) * (p2 ** n2) * (p3 ** n3)
B = 2 
q = 2**10
Bq = 2**9

l = ceil(log(q, B))
sigma = 0.1

# relevant variables
prime, exp = p1, n1
factor = n / (prime**exp)

# rlwe cipher
rlwe = RLWE(n, q, B, sigma=sigma)

# print rlwe-key being used
print("chave: ",rlwe.s)

m = random_msg_gen(B, n)
c = rlwe.enc(m)

base = 2
f = Zx(cyclotomic_polynomial(n))
Zqx = ZZ.quotient(B)['x']
Rb = Zqx.quotient(f)

print("mensagem", m)
for j in Zmod(prime**exp).list_of_elements_of_multiplicative_group():
    i = j - 1
    a , b = c[0], c[1]
    # b = Zx(a.lift())*rlwe.s + m*Bq # calculando em Z[x]

    #print("b calculado", a*rlwe.s + m*Bq)
    #print("b real", rlwe.Rq(b))

    f_new =Zx(Zx(cyclotomic_polynomial(n))(x=x^(i*factor+1)))

    # b_aut = Zx(b(x=x^(i*factor+1))) % Zx(cyclotomic_polynomial(n)) # calculando em Z[x]
    a_aut = Zx(Zx(a.lift())(x=x^(i*factor+1))) % Zx(cyclotomic_polynomial(n))
    b_aut = Zx(Zx(b.lift())(x=x^(i*factor+1)))  % f

    s_aut = Zx(rlwe.s(x=x^(i*factor+1))) % Zx(cyclotomic_polynomial(n))
    
    p = a_aut * s_aut + Zx(m(x=x^(i*factor+1))) * Bq # g(b) = g(a) * g(s) + g(m) q/B + e não está valendo
    print("b aut calculado", rlwe.Rq(p))
    print("b aut real", rlwe.Rq(b_aut))

    base = 2
    K = rlwe.enc_sk(s_aut, base)
    new_cypher = key_switch([a_aut,b_aut], K, base,B, q, n)
    
    print(rlwe.dec(new_cypher))
    m_cmp = Zx(Rb(Zx(m(x=x^(i*factor+1)))).lift())
    print(m_cmp)
    print(rlwe.dec(new_cypher) == m_cmp)

# print("The key switch is working as intended")
