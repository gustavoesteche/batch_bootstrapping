#include <iostream>
#include <NTL/tools.h>

#include "../src/utils/utils.h"
#include "../src/utils/trace.h"
#include "../src/utils/duals.h"
#include "../src/utils/key_switch.h"
#include "../src/utils/RLWE_distribution.h"
#include "../src/schemes/RLWE.h"
#include "../src/schemes/RGSW.h"
#include "../src/packing/cipher_packing.h"

#include "utils_test.h"

using namespace std;

int main() {
    ZZ Q = ZZ(LeftShift(ZZ(1), 100));
    ZZ_p::init(Q);
    
    // validated parameters 
    long p1 = 11, n1 = 1; 
    long p2 = 3, n2 = 2;
    long p3 = 7, n3 = 1;

    //long p1 = 107, n1 = 1; 
    //long p2 = 3, n2 = 1;
    //long p3 = 5, n3 = 1;
    
    long N = pow(p1, n1) * pow(p2, n2) * pow(p3, n3);
    long r = min(euler_phi(pow(p2,n2)), euler_phi(pow(p3,n3)));

    Vec< pair<long, long>> pp; pp.SetLength(3);
    
    pp[0].first = p1; pp[0].second = n1;
    pp[1].first = p2; pp[1].second = n2;
    pp[2].first = p3; pp[2].second = n3;

    ZZX f; ZZ_pXModulus mod_f;
    f = cyclotomic_polynomial(N);
    build(mod_f, conversion_X_to_pX(f));   

    Operations framework = Operations(N, Q, pp, 3.2);

    long rlwe_mode = 1;
    long rgsw_mode = 3;

    Vec<ZZ_pX> MSGS_1 = sample_packages(N, r, pp, mod_f, true);
    pair<Vec<ZZ_pX>, long> p_rlwe = framework.pack_rlwe(MSGS_1, rlwe_mode);

    long tests = 5;
    long depth = 20;
    Vec<ZZ> infty_norm_noises; infty_norm_noises.SetLength(depth/2);
    
    for (int i=0;i<depth/2;i++) infty_norm_noises[i] = ZZ();

    for (long j = 0; j <tests; j++) {
        int k = 0; 
        for (long i = 0; i <depth; i++) { 
            Vec<ZZ_pX> MSGS_2 = sample_packages(N, r, pp, mod_f, true);
            Vec<ZZ_pX> new_MSGS; new_MSGS.SetLength(r);

            for (int i = 0; i < r; i++){
                new_MSGS[i] =  mod_poly((MSGS_1[i] * MSGS_2[i]) % mod_f, 2);
                MSGS_1[i] = new_MSGS[i];
            }

            pair<Mat<ZZ_pX>, long> p_rgsw = framework.pack_rgsw(MSGS_2, rgsw_mode);      
            pair<Vec<ZZ_pX>, long> cipher = framework.ext_prod_trace(p_rlwe, p_rgsw);

            //Vec<ZZ_pX> D_MSGS = framework.unpack_rlwe(cipher);
            if (i%2){
                infty_norm_noises[k] = max(infty_norm_noises[k], framework.get_noise(cipher.first, framework.pack_msgs(new_MSGS, cipher.second)));
                k++;
            }
            
            if (rgsw_mode == 3) {
                rgsw_mode = 4;
            } else {
                rgsw_mode = 3;
            }
            p_rlwe = cipher;

            //cout << D_MSGS << endl;
            //cout << new_MSGS << endl; 
            //cout << compare_vec(D_MSGS, new_MSGS) << endl;  
        }
    }
    cout << infty_norm_noises << endl;
    return 0;
}