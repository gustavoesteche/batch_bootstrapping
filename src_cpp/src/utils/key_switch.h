#ifndef KEYSWITCH_H
#define KEYSWITCH_H


#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include <iostream>

#include "utils.h"

// performs the key-switch in the rlwe ciphertext
Vec<ZZ_pX> ks_rlwe(Vec<ZZ_pX> cipher, const Mat<ZZ_pX>& K, const ZZ& Q, const long l, const long N, const ZZ_pXModulus& mod_phi);

// performs the key-switch in the rgsw ciphertext
Mat<ZZ_pX> ks_gsw(Mat<ZZ_pX> cipher, const Mat<ZZ_pX>& ks_key, const Mat<ZZ_pX>& K, const ZZ& Q, const long l, const long N, const ZZ_pXModulus& mod_phi);

#endif