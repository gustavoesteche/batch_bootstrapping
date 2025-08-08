#ifndef RLWE_H
#define RLWE_H

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <random>
#include <cmath>

#include "../utils/utils.h"
#include "../utils/RLWE_distribution.h"


class RLWE {
public:
    RLWE(const long N, const ZZ& modulus, ZZ_pX sk = ZZ_pX(), double sigma=3.2);

    RLWE() = default;

    ZZ_pX keygen();

    Mat<ZZ_pX> enc_sk(ZZ_pX old_s);

    Vec<ZZ_pX> enc(ZZ_pX msg);

    ZZ_pX dec(Vec<ZZ_pX> cipher);

    ZZ get_noise(Vec<ZZ_pX> cipher, ZZ_pX msg);
    
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
