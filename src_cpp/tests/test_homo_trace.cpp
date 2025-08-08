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
    ZZ Q = ZZ(LeftShift(ZZ(1), 20));
    long l = NumBits(Q-1);

    long p1 = 7, n1 = 1; 
    long p2 = 3, n2 = 2;
    long p3 = 5, n3 = 1;

    long N = pow(p1, n1) * pow(p2, n2) * pow(p3, n3);

    RLWE rlwe = RLWE(N, Q, ZZ_pX(), double(1));
    RGSW rgsw = RGSW(N, Q, rlwe.sk, double(1));

    ZZX f; ZZ_pXModulus mod_f;
    build_poly(f, mod_f, N);
    
    Vec<Vec<long>> autos = construct_aut_power(N, p2, n2);
    Vec<Mat<ZZ_pX>> K_tower = rgsw.gen_evk_k_tensor(p2, n2);
    Mat<ZZ_pX> evk = rgsw.enc(rgsw.sk);

    long tests = 10;
    for (long i = 0; i < tests; i++) {
        ZZ_pX msg = generate_message(N, mod_f);
        Vec<ZZ_pX> cipher = rlwe.enc(msg);

        ZZ_p p_inverse = ZZ_p(1);
        Vec<ZZ_pX> c_trace = homo_generic_trace(cipher, p_inverse, evk, K_tower, Q, l, N, p2, n2, autos, mod_f);
        
        ZZ_pX trace = generic_trace(msg, N, p2, n2, mod_f);
        trace = mod_poly(trace, 2);

        //cout << "msg " << msg << endl;
        //cout << trace << endl;
        //cout << rlwe.dec(c_trace) << endl;
        cout << (rlwe.dec(c_trace) == trace) << endl;
    
    }
    return 0;
}