#include <iostream>

#include "../src/utils/utils.h"
#include "../src/utils/trace.h"
#include "../src/utils/duals.h"
#include "../src/utils/key_switch.h"
#include "../src/utils/RLWE_distribution.h"
#include "../src/schemes/RLWE.h"
#include "../src/schemes/RGSW.h"

#include "utils_test.h"

using namespace std;

int main() {
    ZZ Q = ZZ(LeftShift(ZZ(1), 30));
    long l = NumBits(Q-1);

    long p1 = 31, n1 = 1; 
    long p2 = 3, n2 = 1;
    long p3 = 5, n3 = 1;

    long N = pow(p1, n1) * pow(p2, n2) * pow(p3, n3);

    RLWE rlwe = RLWE(N, Q, ZZ_pX(), double(1.0));
    RGSW rgsw = RGSW(N, Q, rlwe.sk, double(1.0));

    ZZX f; ZZ_pXModulus mod_f;
    build_poly(f, mod_f, N);
    
    long tests = 1;

    for (int i = 0; i < tests; i++) {
        ZZ_pX msg_rlwe = generate_message(N, mod_f);
        ZZ_pX msg_rgsw = generate_message(N, mod_f);

        //cout << msg << endl;
        Vec<ZZ_pX> cipher = rlwe.enc(msg_rlwe);
        Mat<ZZ_pX> evk = rgsw.enc(msg_rgsw);
        //Mat<ZZ_pX> evk = rgsw.enc(mod_poly(rgsw.sk, 2));

        //cout << cipher << endl;
    
        Vec<ZZ_pX> prod_cipher = ext_prod(cipher, evk, Q, l, N, mod_f);
        ZZ_pX res = (msg_rlwe * msg_rgsw) % mod_f;
        
        //sym_mod_poly(res, ZZ(2));
        //cout << mod_poly(res, 2) << endl;
        //cout << rlwe.dec(prod_cipher) << endl;
        //cout << compare_poly(mod_poly(res, 2), rlwe.dec(prod_cipher)) << endl;
        if (!compare_poly(mod_poly(res, 2), rlwe.dec(prod_cipher))) {
            Error("Faulty operation!");
            break;
        }
    }

    cout << "SUCCESS" << endl;
    return 0;
}