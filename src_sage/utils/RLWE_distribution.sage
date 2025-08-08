from sage.stats.distributions.discrete_gaussian_integer \
import DiscreteGaussianDistributionIntegerSampler \
                            as DiscreteGaussian

Zx = PolynomialRing(ZZ, "x")

class RLWE_Distribution:
    def __init__(self, s, N, q, sigma=3.2):
        self.n = N 
        self.f = Zx(cyclotomic_polynomial(N))
        self.s = s 
        self.sigma = sigma
        self.D = DiscreteGaussian(sigma)
        self.Zqx = ZZ.quotient(q)['x'] 
        self.Rq = self.Zqx.quotient(self.f)
    
    def random_noise(self):
        return Zx([self.D() for _ in range(euler_phi(self.n))])
    
    def random_a(self):
        return self.Rq.random_element()
    
    def sample(self):
        s = self.s
        a = self.random_a()
        e = self.random_noise()
        b = (a * s + e) 
        return [a, b]