load("./GSW/GSW.sage")
load("./CRLWE/RLWE.sage")
load("./dual/dbasis.sage")
load("./trace/trace_lib.sage")
load("./trace/rlwe_homo_trace.sage")

Zx = PolynomialRing(ZZ, 'x') 
Qx = PolynomialRing(QQ, 'x')

class Operations:
    def __init__(self, n, q, B, pp, base, sigma):
        self.n, self.q, self.B, self.base = n, q, B, base
        self.l = ceil(log(q, B))
        self.r = min(euler_phi(pp[1][0]**pp[1][1]), euler_phi(pp[2][0]**pp[2][1]))
        self.pp = pp
        self.gsw = GSW(n, q, sigma, B)
        self.rlwe = RLWE(n, q, B, s=self.gsw.sk, sigma=sigma)
        self.m = 1
        for i in range(3):
            self.m *= pp[i][0]**pp[i][1]

        self.Rq = ( ZZ.quotient(q) )['x'].quotient(cyclotomic_polynomial(n))
        self.evk_ktensor, self.auto = [], []
        self.evk = self.gsw.enc(self.gsw.sk)
        for i in range(1,3):
            self.evk_ktensor.append(self.gsw.gen_evk_ktensor_tower(pp[i][0], pp[i][1], base))
            self.auto.append(construct_aut_tower(n, pp[i][0]))

    def trace(self, value, index, exc=0):
        prime, exponent = self.pp[index+1]
        P = prime**exponent
        if exc == 1:    
            P = 1
        res = big_fast_homo_trace_rlwe(value, P, self.evk, self.evk_ktensor[index], self.B, self.q, self.l, self.n, prime, exponent, self.base, self.auto[index])
        return res
         
    def pack_rlwe(self, msgs, mode):
        if mode in [1,2]:
            fator = self.m / (self.pp[mode][0] ** self.pp[mode][1])
            pack = vector(self.rlwe.enc(msgs[0])) * self.Rq(x**(0*fator))
            for i in range(1, self.r): 
                pi = vector(self.rlwe.enc(msgs[i]))
                pack = pack + pi * self.Rq(x**(i*fator))
            return pack, mode

        if mode in [3,4]:
            prime1, exp1 = self.pp[2-mode%2]
            prime2, exp2 = self.pp[1+mode%2]  
            pm = prime1 ** exp1
            db = canon_dbasis(prime1,exp1)
            fator1 = self.m / pm
            fator2 = self.m / (prime2 ** exp2)
            
            pack = vector(self.rlwe.enc(msgs[0])) * Zx(Zx(db[0]*pm)(x=x**fator1)) * self.Rq(x**(0*fator2))
            for i in range(1, self.r): 
                pack = pack + vector(self.rlwe.enc(msgs[i])) * Zx(Zx(db[i]*pm)(x=x**fator1)) * self.Rq(x**(i*fator2))
            return pack, mode

    def unpack_rlwe(self, c):
        msgs = vector(Zx, [0]*self.r)
        mode = c[1]
        cipher = c[0]
    
        if mode in [1,2]:
            prime, exp = self.pp[mode]
            fator = self.m / (prime **exp)
            db = canon_dbasis(prime, exp)
            for i in range(self.r):
                dbi = Zx(Zx(db[i]*(prime**exp))(x = x**fator))    
                ci = self.trace(cipher*dbi, mode-1)
                msgs[i] = self.rlwe.dec(ci)
        
        if mode in [3,4]:
            prime, exp = self.pp[2-mode%2]
            fator = self.m / (prime ** exp)
            index = (mode-1) % 2

            prime1, exp1 = self.pp[1+mode%2]
            fator1 = self.m / (prime1 ** exp1)
            index1 = (mode) % 2
            db = canon_dbasis(prime1, exp1)
            for i in range(self.r):
                ci = self.trace(cipher*Zx(x**(i*fator)), index) 
                dbi = Zx(Zx(db[i]*(prime1**exp1))(x = x**fator1))
                ci = self.trace(ci*dbi,index1)
                msgs[i] = self.rlwe.dec(ci)

        return msgs
 
    def pack_gsw(self, msgs, mode):
        if mode in [1,2]:
            fator = self.m / (self.pp[mode][0] ** self.pp[mode][1])
            pack = Matrix(self.gsw.enc(msgs[0])) * self.Rq(x**(0*fator))
            for i in range(1, self.r): 
                pi = Matrix(self.gsw.enc(msgs[i]))
                pack = pack + pi * self.Rq(x**(i*fator))
            return pack, mode

        if mode in [3,4]:
            prime1, exp1 = self.pp[2-mode%2]
            prime2, exp2 = self.pp[1+mode%2]  
            pm = prime1 ** exp1
            db = canon_dbasis(prime1,exp1)
            fator1 = self.m / pm
            fator2 = self.m / (prime2 ** exp2)
                
            pack = Matrix(self.gsw.enc(msgs[0])) * Zx(Zx(db[0]*pm)(x=x**fator1)) * self.Rq(x**(0*fator2)) 
            for i in range(1, self.r): 
                pack = pack + Matrix(self.gsw.enc(msgs[i])) * Zx(Zx(db[i]*pm)(x=x**fator1)) * self.Rq(x**(i*fator2))
            return pack, mode
    
    def unpack_gsw(self, c):
        raise Exception("Not implemented yet")
        cipher = self.gsw
        return unpack_ciphertext(c[0], self.pp, self.m, self.r, c[1], cipher)
    
    def ext_prod_trace(self, rlwe_c, gsw_c):
        assert rlwe_c[1] != gsw_c[1]
        assert rlwe_c[1] % 2 == gsw_c[1] % 2
        
        index = (rlwe_c[1] - 1) % 2 
        
        result_mode = 2 if rlwe_c[1] % 2 else 1 
        rlwec, gswc = rlwe_c[0], gsw_c[0] 
        exc = rlwe_c[1] // 3
        c = self.gsw.ext_prod(rlwec, gswc)
        return self.trace(c, index, exc), result_mode

    def test_ext_prod(self, rlwe_c, gsw_c):
        return self.gsw.ext_prod(rlwe_c[0], gsw_c[0])

    def add(self, c1, c2):
        assert c1[1] == c2[1], "Different modes"
        return c1[0]+ c2[0], c1[1]

    def mult():
        raise Exception("Not implemented yet.")

# Sample message compatible with multiple product testing
def random_msg_gen_depth(B, n):
    k = ZZ.random_element(0, euler_phi(n))
    k1 = ZZ.random_element(1,2)
    if k == euler_phi(n):
        return 0
    else:
        if k1 == 1:
            return Zx(x^k)
        else:
            return Zx(-x^k)

# Sample the message packages for multiple product testing
def sample_msg_pkgs_depth(B, n,prime_powers, r): 
    M = []
    power = prime_powers[0][0]**prime_powers[0][1]
    cpower = n / power
    for i in range(r):    
        msg = random_msg_gen_depth(B, power)
        msg = Zx(msg(x=x**(cpower)))
        M.append(msg)

    return vector(M)

# Sample random message
def random_msg_gen(B, n):
    f = Zx(cyclotomic_polynomial(n))
    Zbx = ZZ.quotient(B)['x'] 
    Rb = Zbx.quotient(f)
    
    return Zx((Rb.random_element()).lift())

# Sample the message packages
def sample_msg_pkgs(B, n,prime_powers, r): 
    M = []
    power = prime_powers[0][0]**prime_powers[0][1]
    cpower = n / power
    for i in range(r):    
        msg = random_msg_gen(B, power)
        msg = Zx(msg(x=x**(cpower)))
        M.append(msg)

    return vector(M)

# Pack messages according to mode for testing
def pack_msg(packages, prime_powers, m, r, mode, Rq):
    if mode <= 2:
        fator = m / (prime_powers[mode][0] ** prime_powers[mode][1])
        pack = Rq(0)
        for i in range(r): 
            pack = pack + packages[i] * Rq(x**(i*fator))
        return Zx(pack.lift())
    
    db = canon_dbasis(prime_powers[2-mode%2][0],prime_powers[2-mode%2][1])
    pm = prime_powers[2-mode%2][0]**prime_powers[2-mode%2][1]
    fator1 = m / (prime_powers[2-mode%2][0] ** prime_powers[2-mode%2][1])
    fator2 = m / (prime_powers[1+mode%2][0] ** prime_powers[1+mode%2][1])
        
    pack = Rq(0) 
    for i in range(r): 
        pack = pack + packages[i] * Rq(x**(i*fator2)) * Zx(Zx(db[i]*pm)(x=x**fator1))
    return Zx(pack.lift())