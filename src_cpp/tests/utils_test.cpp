#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <iostream>
#include <random>

using namespace NTL;
using namespace std;

#include "../src/utils/utils.h"


long get_random_long(long upper_limit) {
    std::random_device rd;
    std::mt19937_64 gen(rd()); // 64-bit Mersenne Twister
    std::uniform_int_distribution<long> dis(0, upper_limit);
    return dis(gen);
}


void build_poly(ZZX &f, ZZ_pXModulus &mod_f, const long N) { 
    f = cyclotomic_polynomial(N);
    build(mod_f, conversion_X_to_pX(f));   
}

ZZ_pX generate_message(long N, const ZZ_pXModulus& mod_phi) {
    ZZ_pX msg; msg.rep.SetLength(N);
    long phi_n = euler_phi(N);
    for (long i = 0; i < phi_n; ++i) {
        msg.rep[i] = (conv<ZZ>(random_ZZ_p()) % 2); 
    }
    rem(msg, msg, mod_phi);
    return msg;   
}

ZZ_pX generate_message_depth(long N) {
    long phi_n = euler_phi(N);    
    ZZ_pX msg; msg.rep.SetLength(N);
    SetCoeff(msg, get_random_long(phi_n-1), 1);
    return msg;
}

Vec<ZZ_pX> sample_packages(long N, long r, Vec<pair<long, long>> pp, const ZZ_pXModulus& mod_phi, bool depth=false) {
    Vec<ZZ_pX> M; M.SetLength(r);
    long power = pow(pp[0].first, pp[0].second); long cpower = N / power;
    if (depth) {
        for (int i=0; i<r; i++)
            M[i] = aut(generate_message_depth(power), cpower, N, mod_phi);
        
        return M;
    }

    for (int i=0; i<r; i++)
        M[i] = aut(generate_message(power, mod_phi), cpower, N, mod_phi);
    return M;
}

bool compare_poly(const ZZ_pX& p1, const ZZ_pX& p2) {
    if (deg(p1) != deg(p2)) 
        return false;

    for (int i = 0; i < deg(p1); i++){
        if (p1[i] != p2[i])
            return false;
    }
    return true;
}

bool compare_poly(const ZZX& p1, const ZZX& p2) {
    if (deg(p1) != deg(p2)) 
        return false;

    for (int i = 0; i < deg(p1); i++){
        if (p1[i] != p2[i])
            return false;
    }
    return true;
}

bool compare_vec(const Vec<ZZ_pX>& v1, const Vec<ZZ_pX>& v2) { 
    if (v1.length() != v2.length()) 
        return false; 

    for (long i = 0; i < v1.length() ; i++) {
        if (!compare_poly(v1[i], v2[i]))
            return false;
    }
    
    return true;
}