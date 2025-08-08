#include <iostream>

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
    ZZ Q = ZZ(LeftShift(ZZ(1), 20));
    long l = NumBits(Q-1);
    ZZ_p::init(Q);
    long p1 = 31, n1 = 1; 
    long p2 = 3, n2 = 1;
    long p3 = 5, n3 = 1;

    long N = pow(p1, n1) * pow(p2, n2) * pow(p3, n3);
    long r = min(euler_phi(pow(p2,n2)), euler_phi(pow(p3,n3)));

    Vec< pair<long, long>> pp; pp.SetLength(3);
    
    pp[0].first = p1; pp[0].second = n1;
    pp[1].first = p2; pp[1].second = n2;
    pp[2].first = p3; pp[2].second = n3;
    
    Operations framework = Operations(N, Q, pp, 1.0);

    ZZX f; ZZ_pXModulus mod_f;
    f = cyclotomic_polynomial(N);
    build(mod_f, conversion_X_to_pX(f));   

    //ZZ_pX msg = generate_message(N, mod_f);
    //Vec<ZZ_pX> c = framework.rlwe.enc(msg);
    //cout << mod_poly(generic_trace(msg, N, p2, n2, mod_f), 2) << endl;
    //cout << framework.rlwe.dec(framework.trace(c, 0, true)) << endl;
    long tests = 100;
    for (long i = 0; i <tests; i++) { 
        Vec<ZZ_pX> MSGS = sample_packages(N, r, pp, mod_f);
        pair<Vec<ZZ_pX>, long> p_rlwe = framework.pack_rlwe(MSGS, 1);
        Vec<ZZ_pX> D_MSGS = framework.unpack_rlwe(p_rlwe);
        Vec<ZZ_pX> MSGS_; MSGS_.SetLength(r);
        
        for (int i = 0; i<r; i++) {
            MSGS_[i] = mod_poly(MSGS[i], 2);
        } 
        
        //cout << MSGS_ << endl;
        //cout << D_MSGS << endl;
        cout << compare_vec(MSGS_, D_MSGS) << endl;
    }
    return 0;
}