#include "utils.h"
#include <NTL/ZZ_pXFactoring.h>
#include <cmath>

#include "RLWE_distribution.h"

using namespace std;
using namespace NTL;

DiscreteGaussian::DiscreteGaussian(double sigma)
    : distribution(0.0, sigma) {}

ZZ DiscreteGaussian::sample() {
    return ZZ(lround(distribution(generator)));
}

RLWE_Distribution::RLWE_Distribution(const ZZ_pX& sk, long N, const ZZ& modulus, double sigma)
    : n(N), q(modulus), s(sk), D(sigma)
{
    ZZ_p::init(q);
    phi_n = euler_phi(n);
    f = cyclotomic_polynomial(N);
    ZZ_pX f_p = conversion_X_to_pX(f);
    build(mod_f, f_p);
}

ZZ_pX RLWE_Distribution::random_noise() {
    ZZ_pX e;
    e.SetMaxLength(phi_n-1);
    for (long i = 0; i <= phi_n; ++i) {
        SetCoeff(e, i, to_ZZ_p(D.sample()));
    }
    return e;
}

ZZ_pX RLWE_Distribution::random_a() {
    ZZ_pX a;
    random(a, deg(f));
    rem(a, a, mod_f);
    return a;
}

Vec<ZZ_pX> RLWE_Distribution::sample() {
    ZZ_pX a = random_a();
    ZZ_pX e = random_noise();
    ZZ_pX b;
    
    b = (a * s) % mod_f;
    add(b, b, e);
    rem(b, b, mod_f);

    Vec<ZZ_pX> ct; ct.SetLength(2);
    ct[0] = a; ct[1] = b;
    return ct;
}
