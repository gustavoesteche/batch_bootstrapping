#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include <iostream>

#include "utils.h"

using namespace std;
using namespace NTL;

Vec<ZZ_pX> ks_rlwe(Vec<ZZ_pX> cipher, const Mat<ZZ_pX>& K, const ZZ& Q, const long l, const long N, const ZZ_pXModulus& mod_phi) {
    Vec<ZZ_pX> evak;
    Vec<ZZ_pX> v; v.SetLength(2);
    v[0] = ZZ_pX(); v[1] = cipher[1];

    Mul(evak, inv_g_poly_p(cipher[0], Q, l, N), K, mod_phi);
    v[0] = - evak[0]; v[1] = v[1] - evak[1];

    return v;
}

Mat<ZZ_pX> ks_gsw(Mat<ZZ_pX> cipher, const Mat<ZZ_pX>& ks_key, const Mat<ZZ_pX>& K, const ZZ& Q, const long l, const long N, const ZZ_pXModulus& mod_phi) {
    Mat<ZZ_pX> new_cipher;
    new_cipher.SetDims(2*l, 2);

    for (int i = 0; i<l; i++) {
        new_cipher[l+i] = ks_rlwe(cipher[l+i], K, Q, l, N, mod_phi);
        Mul(new_cipher[i],inv_g_ciphertext_p(new_cipher[l+i], Q, l, N), ks_key, mod_phi);
    }

    return new_cipher;
}
