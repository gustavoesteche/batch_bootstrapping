#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <cmath>
#include <NTL/RR.h>

#include "../utils/utils.h"
#include "../utils/RLWE_distribution.h"
#include "RGSW.h"

using namespace std;
using namespace NTL;

RGSW::RGSW(const long N, const ZZ& modulus, ZZ_pX sk, double sigma)
    : N(N), Q(modulus), sk(sk), sigma(sigma)
{
    l = NumBits(Q-1);
    Q = power(ZZ(2), l);
    ZZ_p::init(Q);
    
    phi_N = euler_phi(N);
    f = cyclotomic_polynomial(N);
    ZZ_pX f_p = conversion_X_to_pX(f);
    build(mod_f, f_p);
    
    if (sk == ZZ_pX())
        sk = keygen();
    D = RLWE_Distribution(sk, N, Q, sigma);
}

ZZ_pX RGSW::keygen() {
    while (sk == ZZ_pX()) {
        random(sk, phi_N-1);
        rem(sk, sk, mod_f);
    }
    return sk;
}

void RGSW::enc_sk(Mat<ZZ_pX> & mat_c, ZZ_pX old_s){
    ZZ_p power_up = ZZ_p(1);
    for (int i = 0; i < l; i++) { 
        mat_c[i] = D.sample();
        mat_c[i][1] = mat_c[i][1] + old_s * power_up;   
        power_up = power_up * 2;
    }
}

Vec<Mat<ZZ_pX>> RGSW::gen_evk_k_tensor(const long prime, const long exp){
    long size_ = euler_phi(conv<long>(power(ZZ(prime), exp)));
    long m = N;
    Vec<Mat<ZZ_pX>> res; res.SetLength(size_);
    long k = 0;
    for (long i = 0; i<exp; i++){
        Vec<long> autos = find_aut_tower(m, prime);
        for (auto x: autos) {
            res[k].SetDims(l, 2);
            enc_sk(res[k], aut(sk, x, N, mod_f));          
            k++;
        }

        m /= prime;
    }
    return res;
}

Mat<ZZ_pX> RGSW::enc(ZZ_pX msg) {
    Mat<ZZ_pX> C; C.SetDims(2*l, 2);
    rem(msg, msg, mod_f);

    //msg = mod_poly(msg, 2);
    ZZ_p power_up = ZZ_p(1);
    for (int j = 0; j < 2;j++) {
        for (int i = j*l; i<(j+1)*l; i++) {
            C[i] = D.sample();
            C[i][j] = C[i][j] + msg * power_up; 
            power_up *= 2;
        }
        power_up = ZZ_p(1);
    }
    
    return C;
}

ZZ_pX RGSW::dec(Mat<ZZ_pX> cipher) {
    ZZ_pX noisy_msg = cipher[2*l-1][1] - (cipher[2*l-1][0] * sk);
    rem(noisy_msg, noisy_msg, mod_f);

    ZZ_pX msg; clear(msg);
    RR scale = to_RR(2) / to_RR(Q);  

    for (long i = 0; i <= deg(f); ++i) {
        ZZ_p c = (i <= deg(noisy_msg)) ? noisy_msg[i] : ZZ_p(0);
        RR val = to_RR(conv<ZZ>(c)) * scale;
        ZZ rounded = RoundToZZ(val);
        sym_mod(rounded, ZZ(2));
        SetCoeff(msg, i, to_ZZ_p(rounded));
    }
    return msg;
}

Vec<ZZ_pX> RGSW::inv_g_row_ciphertext(Vec<ZZ_pX> c) {
    Vec<ZZ_pX> res; res.SetLength(2*l);
    for (int i = 0; i < 2; i++) {
        Vec<ZZ_pX> ans = inv_g_poly_p(c[i], Q, l, N);
        for (int j = 0; j<l; j++) {
            res[i*l+j] = ans[j];
        }
    }
    return res;
}

Mat<ZZ_pX> RGSW::add(Mat<ZZ_pX> c1, Mat<ZZ_pX> c2) {
    Mat<ZZ_pX> res; res.SetDims(2*l,2);
    for (int i=0;i<2*l;i++) { 
        for (int j=0;j<2;j++){
            clear(res[i][j]);
            res[i][j] = c1[i][j] + c2[i][j];
            rem(res[i][j], res[i][j], mod_f);
        }
    }
    return res;
}

Mat<ZZ_pX> RGSW::mult(Mat<ZZ_pX> c1, Mat<ZZ_pX> c2) {
    Mat<ZZ_pX> res; res.SetDims(2*l, 2);
    for (int i = 0; i<2*l; i++) {
        Mul(res[i], inv_g_row_ciphertext(c1[i]), c2, mod_f);
    }
    return res;
}

Vec<ZZ_pX> RGSW::ext_prod(Vec<ZZ_pX> rlwe_c, Mat<ZZ_pX> gsw_c) {
    Vec<ZZ_pX> res; res.SetLength(2);
    Mul(res, inv_g_row_ciphertext(rlwe_c), gsw_c, mod_f);
    return res;
}

ZZ RGSW::get_noise(Mat<ZZ_pX> cipher, ZZ_pX msg) {
    
}