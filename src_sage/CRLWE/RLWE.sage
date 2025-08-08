from sage.stats.distributions.discrete_gaussian_integer \
import DiscreteGaussianDistributionIntegerSampler \
                            as DiscreteGaussian

load("./utils/utils.sage")
load("./utils/RLWE_distribution.sage")
load("./utils/gadget_matrix.sage")

# RLWE with cipher and cipher functioning like a GSW row 
class RLWE:
    def __init__(self, N, q, B, s = None, sigma=3.2):
        self.n = N 
        self.f = cyclotomic_polynomial(N)
        self.sigma = sigma
        self.B = B
        self.q = q 
        self.l = ceil(log(q, B))
        self.Zqx = ZZ.quotient(q)['x'] 
        self.Rq = self.Zqx.quotient(self.f)
        self.Rq = ( ZZ.quotient(q) )['x'].quotient(self.f)
        self.Rb = ( ZZ.quotient(B) )['x'].quotient(self.f)
        if s == None:
            s = self.keygen()
        self.s = s
        self.dist_rlwe = RLWE_Distribution(self.s, N, q, sigma)
    
    def keygen(self):
        sk = 0
        while 0 == sk:
            sk = Zx([ZZ.random_element(0, self.q) for _ in range(euler_phi(self.n))])
        return sk

    def enc_sk(self, sk, base):
        c_sk = Matrix(Zx, self.l, 2)
        for i in range(self.l):
            c = self.dist_rlwe.sample()
            c_sk[i] = [Zx(c[0].lift()), Zx(c[1].lift()) + sk*(base**i)]

        return c_sk

    def enc(self, msg):
        c = self.dist_rlwe.sample()
        factor = (self.q / self.B)
        return [c[0], c[1] + self.Zqx(msg * factor)] 

    def dec(self, ciphertext):
        noisy_m = Zx((ciphertext[1] - self.s * ciphertext[0]).lift())
        rounded_m = round_poly(noisy_m * self.B/self.q)
        return sym_mod_poly(rounded_m , B)

    def get_noise(self, ciphertext, msg=None):
        if msg == None:
            msg = self.dec(ciphertext)
        
        factor = (self.q / self.B)
        error = Zx((ciphertext[1] - ciphertext[0]*self.s - self.Zqx(msg * factor)).lift())
        error = sym_mod_poly(error, self.q)

        return infinity_norm_vec(error) 

    def get_log_noise(self, ciphertext, msg=None):
        norm_e = self.get_noise(ciphertext, msg)
        if norm_e == 0:
            return -1
        return log(norm_e, 2).n() 
        
def random_msg_gen(B, n):
    return Zx(sum([ZZ.random_element(0, B) * x**i for i in range(euler_phi(n)-1)])) % Zx(cyclotomic_polynomial(n))
    