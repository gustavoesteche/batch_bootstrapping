#include <iostream>

#include "../src/utils/utils.h"
#include "../src/utils/trace.h"
#include "../src/utils/duals.h"
#include "../src/utils/key_switch.h"
#include "../src/utils/RLWE_distribution.h"
#include "../src/schemes/RLWE.h"

#include "utils_test.h"

using namespace std;

int main() {
    ZZ Q = ZZ(LeftShift(ZZ(1), 20));
    long l = NumBits(Q-1);

    long p1 = 1, n1 = 1; 
    long p2 = 3, n2 = 1;
    long p3 = 5, n3 = 1;

    long N = pow(p1, n1) * pow(p2, n2) * pow(p3, n3);

    RLWE old_rlwe = RLWE(N, Q, ZZ_pX(), double(2));
    RLWE new_rlwe = RLWE(N, Q, ZZ_pX(), double(2));
    
    ZZX f; ZZ_pXModulus mod_f;
    build_poly(f, mod_f, N);
    
    ZZ_pX msg = generate_message(N, mod_f);
    
    Vec<ZZ_pX> cipher = old_rlwe.enc(msg);
    
    ZZ_pX tst = old_rlwe.dec(cipher);
    Mat<ZZ_pX> K = new_rlwe.enc_sk(old_rlwe.sk);
    Vec<ZZ_pX> new_cipher = ks_rlwe(cipher, K, Q, l, N, mod_f);

    cout << msg << endl;
    cout << compare_poly(new_rlwe.dec(new_cipher), mod_poly(msg, 2)) << endl;
    return 0;
}