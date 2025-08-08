#ifndef RGSW_H
#define RGSW_H

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <random>
#include <cmath>

#include "../utils/utils.h"
#include "../utils/RLWE_distribution.h"

class RGSW {
public:
    RGSW(const long N, const ZZ& modulus, ZZ_pX sk = ZZ_pX(), double sigma=3.2);
    RGSW() = default;

    ZZ_pX keygen();

    void enc_sk(Mat<ZZ_pX>& mat_c, ZZ_pX old_s);

    Vec<Mat<ZZ_pX>> gen_evk_k_tensor(const long prime, const long exp);

    Mat<ZZ_pX> enc(ZZ_pX msg);

    ZZ_pX dec(Mat<ZZ_pX> cipher);

    Vec<ZZ_pX> inv_g_row_ciphertext(Vec<ZZ_pX> c);
    
    Mat<ZZ_pX> add(Mat<ZZ_pX> c1, Mat<ZZ_pX> c2);

    Mat<ZZ_pX> mult(Mat<ZZ_pX> c1, Mat<ZZ_pX> c2);

    Vec<ZZ_pX> ext_prod(Vec<ZZ_pX> rlwe_c, Mat<ZZ_pX> gsw_c);

    ZZ get_noise(Mat<ZZ_pX> cipher, ZZ_pX msg);
    
    ZZ_pX sk;
    
private:
    long N;
    long phi_N;
    ZZ Q;
    long l;
    double sigma;
    
    ZZX f;
    ZZ_pXModulus mod_f;
    RLWE_Distribution D;
};

#endif
