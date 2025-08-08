#include <iostream>

#include "../src/utils/utils.h"
#include "../src/utils/trace.h"
#include "../src/utils/duals.h"
#include "../src/utils/key_switch.h"
#include "../src/utils/RLWE_distribution.h"
#include "../src/schemes/RLWE.h"
#include "../src/schemes/RGSW.h"

using namespace std;

int main() {
    ZZ Q = ZZ(7001);
    long l = NumBits(Q - 1) + 1;
    long p = 7, n = 1; long N=1;
    for (int i = 0; i<n; i++)
        N *= p;

    RGSW rgsw = RGSW(N, Q, ZZ_pX(), double(2));
    
    ZZ_pX msg;
    long deg = 6;
    msg.rep.SetLength(deg);
    for (long i = 0; i < deg; ++i) {
        msg.rep[i] = to_ZZ_p(conv<ZZ>(random_ZZ_p()) % 2); 
    }
    msg.normalize();
    
    ZZX f = cyclotomic_polynomial(N);
    ZZ_pXModulus mod_f;
    ZZ_pX f_p = conversion_X_to_pX(f);
    build(mod_f, f_p);
    rem(msg, msg, mod_f);
    
    Mat<ZZ_pX> cipher = rgsw.enc(msg);
    ZZ_pX tst = rgsw.dec(cipher);

    cout << msg << endl;
    cout << tst << endl;
    return 0;
}