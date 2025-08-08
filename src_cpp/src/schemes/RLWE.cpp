#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <cmath>
#include <NTL/RR.h>

#include "../utils/utils.h"
#include "../utils/RLWE_distribution.h"
#include "RLWE.h"

using namespace std;
using namespace NTL;

RLWE::RLWE(const long N, const ZZ& modulus, ZZ_pX sk, double sigma)
    : N(N), Q(modulus), sk(sk), sigma(sigma)
{
    l = NumBits(Q-1);
    Q = power(ZZ(2), l);
    ZZ_p::init(Q);

    // aqui o problema
    phi_N = euler_phi(N);
    f = cyclotomic_polynomial(N);
    ZZ_pX f_p = conversion_X_to_pX(f);
    build(mod_f, f_p);
    
    if (sk == ZZ_pX())
        sk = keygen();
    D = RLWE_Distribution(sk, N, Q, sigma);
}

ZZ_pX RLWE::keygen() {
    while (sk == ZZ_pX()) {
        random(sk, phi_N-1);
        rem(sk, sk, mod_f);
    }
    return sk;
}

Mat<ZZ_pX> RLWE::enc_sk(ZZ_pX old_s) {
    Mat<ZZ_pX> mat_c; mat_c.SetDims(l, 2); 
    Vec<ZZ_pX> c;
    ZZ_p power_up = ZZ_p(1);
    for (int i = 0; i < l; i++) { 
        c = D.sample();
        c[1] = c[1] + old_s * power_up;
        mat_c[i] = c;   
        power_up = power_up * 2;
    }
    return mat_c;
}

Vec<ZZ_pX> RLWE::enc(ZZ_pX msg) {
    Vec<ZZ_pX> cipher = D.sample();
    msg = mod_poly(msg, 2);
    cipher[1] = cipher[1] + msg * to_ZZ_p(Q/2);
    return cipher;
}

ZZ_pX RLWE::dec(Vec<ZZ_pX> cipher) {
    ZZ_pX noisy = cipher[1] - (cipher[0] * sk);
    rem(noisy, noisy, mod_f);

    ZZ_pX msg;
    RR factor = to_RR(2) / to_RR(Q);

    for (long i = 0; i <= deg(f); ++i) {
        ZZ coeff_z = (i <= deg(noisy)) ? rep(noisy[i]) : ZZ(0);
        RR scaled = to_RR(coeff_z) * factor;
        ZZ rounded = RoundToZZ(scaled);
        sym_mod(rounded, ZZ(2));
        SetCoeff(msg, i, to_ZZ_p(rounded));
    }
    return msg;
}

ZZ RLWE::get_noise(Vec<ZZ_pX> cipher, ZZ_pX msg) { 
    ZZ_pX noise = cipher[1] - (cipher[0] * sk) % mod_f - msg * to_ZZ_p(Q/2);

    return infinity_norm(sym_mod_poly(noise, Q));
}
