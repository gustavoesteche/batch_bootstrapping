#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <iostream>

#include "utils.h"
#include "key_switch.h"

using namespace std;
using namespace NTL;


ZZ_pX generic_trace(const ZZ_pX& poly, long m, long  p, long n, const ZZ_pXModulus& mod_phi) {
    Vec<long> auts = find_aut(m, p, n);
    ZZ_pX ans = ZZ_pX(0);
    for (int i = 0; i < auts.length(); i++) {
        ans = ans + aut(poly, auts[i], m, mod_phi); 
    }

    return ans;
}

Vec<ZZ_pX> homo_generic_trace(const Vec<ZZ_pX>& rlwe_cipher, ZZ_p p_inverse, Mat<ZZ_pX> evk, Vec<Mat<ZZ_pX>> K_tower, 
    const ZZ Q, const long l, const long N, const long prime, const long exp, Vec< Vec<long> > autos, const ZZ_pXModulus& mod_phi){
    
    Vec<ZZ_pX> rlwe_c; rlwe_c.SetLength(2);
    rlwe_c[0] = ZZ_pX(); rlwe_c[1] = rlwe_cipher[0] * p_inverse;
    Vec<ZZ_pX> d = ext_prod(rlwe_c, evk, Q, l, N, mod_phi);
    
    long k = 0;
    Vec<ZZ_pX> ks_res; ks_res.SetLength(2);
    for (int i = 0; i<exp;i++){
        clear(rlwe_c[0]); clear(rlwe_c[1]);
        for (int j = 0; j < autos[i].length();j++) {
            ks_res[0] = aut(d[0], autos[i][j], N, mod_phi);
            ks_res[1] = aut(d[1], autos[i][j], N, mod_phi);
            ks_res = ks_rlwe(ks_res, K_tower[k], Q, l, N, mod_phi);
            rlwe_c[0] += ks_res[0];rlwe_c[1] += ks_res[1];
            k++;
        }
        d = rlwe_c;
    }
    ks_res[0] = - d[0];
    ks_res[1] = generic_trace(rlwe_cipher[1], N, prime, exp, mod_phi) * p_inverse - d[1];
    
    return ks_res;
}
