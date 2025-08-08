load("./utils/utils.sage")
load("./utils/RLWE_distribution.sage")
load("./utils/gadget_matrix.sage")

Zx = ZZ['x']

class GSW:
    def __init__(self, n, q, sigma, B = 2):
        self.n, self.sigma, self.B = n, sigma, B
        self.f = Zx(cyclotomic_polynomial(n))
        self.l = ceil(log(q, B))
        self.q = B^self.l
        g = Matrix(ZZ, self.l, 1, [B^i for i in range(self.l)])
        I = Matrix.identity(2)
        self.G = I.tensor_product(g) 
        self.Rq = ( ZZ.quotient(q) )['x'].quotient(self.f)
        self.Rb = ( ZZ.quotient(B) )['x'].quotient(self.f)
        self.sk = self.keygen()
        self.dist_rlwe = RLWE_Distribution(self.sk, n, q, sigma)
        self.Ks = self.enc(-self.sk)
    
    def keygen(self):
        sk = 0
        while 0 == sk:
            sk = Zx([ZZ.random_element(0, self.q) for _ in range(euler_phi(self.n))])
        return sk
    
    def enc(self, m):
        Rq, l, sk = self.Rq, self.l, self.sk
        C = Matrix(Rq, 2* l, 2)
        for i in range(2* l):
            C[i] = self.dist_rlwe.sample()
     
        C += m * self.G 
        return C
    
    def dec(self, C):
        B, l, q, sk = self.B, self.l, self.q, self.sk
        a = C[2* l - 1, 0]
        b = C[2* l - 1, 1] 
        noisy_m = Zx((b - a* sk).lift())
        rounded_m = round_poly(B * noisy_m / q) 
        return sym_mod_poly(rounded_m , B)

    def enc_sk(self, sk, base):
        c_sk = Matrix(Zx, self.l, 2)
        for i in range(self.l):    
            c = self.dist_rlwe.sample()
            c_sk[i] = [Zx(c[0].lift()), Zx(c[1].lift()) + (sk)*(base**i)]

        return c_sk

    def gen_evk_ktensor(self, prime, exp, base):
        mi = prime**exp
        factor = self.n / mi
        evk_ktensor = []
        auto = find_aut(self.n, prime, exp)

        for i in auto:
            sk_t = big_aut(self.sk, i, Zx) % self.f
            evk_ktensor.append(self.enc_sk(sk_t, base))

        return evk_ktensor

    def gen_evk_ktensor_tower(self, prime, exp, base):
        evk_ktensor = []
        m = self.n
        for j in range(exp):
            auto = find_aut_tower(m, prime)
            K = []
            for i in auto:
                sk_t = big_aut(self.sk, i, Zx) % self.f
                K.append(self.enc_sk(sk_t, base))
            m = m / prime
            evk_ktensor.append(K)
        return evk_ktensor
            
    def add(self, C0, C1):
        return  C0 + C1 

    def inv_g_row_ciphertext(self, c):
        B, l, n = self.B, self.l, self.n
        a, b = c[0], c[1]
        res = vector(Zx, [0] * 2 * l)
        res[0:l] = inv_g_poly(Zx(a.lift()), B, self.q, n)
        res[l:2 * l] = inv_g_poly(Zx(b.lift()), B, self.q, n)
        return res
 
    def mult(self, C0, C1):
        result = Matrix(self.Rq, 2 * self.l, 2)

        for i in range(2 * self.l):
            decomp = self.inv_g_row_ciphertext(C0[i])
            prod_in_Rq = decomp * C1
            result[i] = prod_in_Rq
        return result

    def ext_prod(self, rlwe_c, gsw_c):
        return self.inv_g_row_ciphertext(rlwe_c) * gsw_c

    def get_noise(self, C, msg=None):
        l, sk = self.l, self.sk
        if msg == None:
            msg = self.dec(C) 
        
        C -= msg * self.G
        a = vector(C[:, 0])
        b = vector(C[:, 1])
        e = b - a*sk
        e = sym_mod_vec(e, self.q) 
        
        return infinity_norm_vec(e)
    
    def get_log_noise(self, C, msg=None):
        norm_e = self.get_noise(C, msg)
        if norm_e == 0:
            return -1
        return log(norm_e, 2).n() 
