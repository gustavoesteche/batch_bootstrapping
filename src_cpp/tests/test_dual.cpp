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
    ZZ Q = ZZ(LeftShift(ZZ(1), 10));
    long l = NumBits(Q-1);
    ZZ_p::init(Q);
    long p1 = 3, n1 = 2; 
    long p2 = 2, n2 = 1;
    long p3 = 5, n3 = 1;

    long N = pow(p1, n1) * pow(p2, n2) * pow(p3, n3);

    ZZX f; ZZ_pXModulus mod_f;
    build_poly(f, mod_f, N);
    
    // mensagem
    ZZ_pX poly; 
    ZZ power = conv<ZZ>(pow(p1, n1));
    long fator = N / pow(p1, n1); 

    Vec<ZZ_pX> db = canon_dbasis(p1, n1, mod_f);
    ZZ_pX dbfator;
    for (int i=0;i<db.length();i++) {
        clear(poly); SetCoeff(poly, i*fator, 1);
        dbfator = aut(db[i], fator, N, mod_f);
        ZZ_pX a = poly * dbfator * to_ZZ_p(InvMod(power, Q)) % mod_f; 
        cout << generic_trace(a, N, p1, n1, mod_f) << endl;
    }
    return 0;
}