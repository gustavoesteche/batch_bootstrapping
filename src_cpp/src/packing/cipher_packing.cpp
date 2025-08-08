#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include <NTL/RR.h>
#include <iostream>

#include "../utils/utils.h"
#include "../utils/trace.h"
#include "../utils/duals.h"
#include "../utils/key_switch.h"
#include "../utils/RLWE_distribution.h"
#include "../schemes/RLWE.h"
#include "../schemes/RGSW.h"

#include "cipher_packing.h"


Operations::Operations(const long &N, const ZZ& modulus, Vec<pair<long, long>> pp, double sigma)
    : N(N), Q(modulus), pp(pp)
{
    l = NumBits(Q-1);
    Q = power(ZZ(2), l);
    ZZ_p::init(Q);

    powers.push_back(pow(pp[1].first, pp[1].second));
    powers.push_back(pow(pp[2].first, pp[2].second));
    r = min(euler_phi(powers[0]), euler_phi(powers[1]));
    phi_N = euler_phi(N);
    f = cyclotomic_polynomial(N);
    build(mod_f, conv<ZZ_pX>(f));
    rlwe = RLWE(N, Q, ZZ_pX(), sigma);
    rgsw = RGSW(N, Q, rlwe.sk, sigma);
    evk = rgsw.enc(rgsw.sk);
    gen_eval_keys();
    gen_canon_dbasis();
    construct_aut();
}

void Operations::gen_canon_dbasis() {
    for (int i = 1; i<=2; i++) {
        vec_dbasis.push_back(canon_dbasis(pp[i].first, pp[i].second, mod_f));
    }    
}

void Operations::gen_eval_keys() {
    evk_ktensor.clear();
    for (long i = 1; i <= 2; ++i) {
        auto kts = rgsw.gen_evk_k_tensor(pp[i].first, pp[i].second);
        evk_ktensor.push_back(kts);
    }
}       

void Operations::construct_aut() {
    auto_ks.clear();
    for (long i = 1; i <= 2; ++i) {
        auto aut = construct_aut_power(N, pp[i].first, pp[i].second);
        auto_ks.push_back(aut);
    }
}

Vec<ZZ_pX> Operations::trace(Vec<ZZ_pX> poly, long index) {
    long prime = pp[index + 1].first, exp = pp[index + 1].second;    
    return homo_generic_trace(poly, ZZ_p(1), evk, evk_ktensor[index], Q, l, N, prime, exp, auto_ks[index], mod_f);
}

pair<Vec<ZZ_pX>, long> Operations::ext_prod_trace(pair<Vec<ZZ_pX>, long> p_rlwe, pair<Mat<ZZ_pX>, long> p_rgsw) {
    if (p_rlwe.second == p_rgsw.second || p_rgsw.second % 2 != p_rlwe.second % 2) {
        Error("Operands dont match");
    }

    long index = (p_rlwe.second - 1) % 2;
    long result_mode = p_rlwe.second % 2 == 1 ? 2 : 1;
    Vec<ZZ_pX> mult = rgsw.ext_prod(p_rlwe.first, p_rgsw.first);
    mult[0] = mult[0] * to_ZZ_p(InvMod(conv<ZZ>(powers[index]), Q));
    mult[1] = mult[1] * to_ZZ_p(InvMod(conv<ZZ>(powers[index]), Q));
    return {trace(mult, index), result_mode};
}

pair<Vec<ZZ_pX>, long> Operations::pack_rlwe(const Vec<ZZ_pX>& msgs, const long& mode) {
    Vec<ZZ_pX> pack; pack.SetLength(2); clear(pack[0]); clear(pack[1]);
    
    if (mode == 1 || mode == 2) {
        long fator = N / powers[mode-1];
        Vec<ZZ_pX> sum; ZZ_pX power;

        for (long i = 0; i < r; i++) {
            clear(power); SetCoeff(power, i*fator, 1);
            sum = rlwe.enc(msgs[i]);
            pack[0] = pack[0] + (sum[0] * power) % mod_f;
            pack[1] = pack[1] + (sum[1] * power) % mod_f;
        }
        return {pack, mode};
    }

    if (mode == 3 || mode == 4) {
        long fator1 = N / powers[1 - mode % 2];
        long fator2 = N / powers[mode % 2];

        Vec<ZZ_pX> sum; ZZ_pX power; ZZ_pX db_fator;
        for (long i = 0; i < r; ++i) {
            clear(power); SetCoeff(power, i*fator2, 1);
            db_fator = aut(vec_dbasis[mode - 3][i], fator1, N, mod_f);
            sum = rlwe.enc(msgs[i]);
            pack[0] = pack[0] + (((sum[0] * power) % mod_f) * db_fator) % mod_f;
            pack[1] = pack[1] + (((sum[1] * power) % mod_f) * db_fator) % mod_f;
        }
        return {pack, mode};
    }

    Error("Operation mode not implemented");
}

Vec<ZZ_pX> Operations::unpack_rlwe(const pair<Vec<ZZ_pX>, long>& p_rlwe) {
    Vec<ZZ_pX> MSGS; MSGS.SetLength(r);
    
    if (p_rlwe.second == 1 || p_rlwe.second == 2) {
        long fator = N / powers[p_rlwe.second-1];
        
        ZZ_pX dbfator;
        Vec<ZZ_pX> c; c.SetLength(2);
        for (int i = 0; i<r; i++) {
            dbfator = aut(vec_dbasis[p_rlwe.second-1][i], fator, N, mod_f);

            c[0] = (p_rlwe.first[0] * dbfator) % mod_f * to_ZZ_p(InvMod(conv<ZZ>(powers[p_rlwe.second-1]), Q)); 
            c[1] = (p_rlwe.first[1] * dbfator) % mod_f * to_ZZ_p(InvMod(conv<ZZ>(powers[p_rlwe.second-1]), Q));
            c = trace(c, p_rlwe.second-1);
            MSGS[i] = rlwe.dec(c);
        }
        return MSGS;
    }

    if (p_rlwe.second == 3 || p_rlwe.second == 4) {
        long fator1, fator2;
        
    } 

    Error("Operation mode not implemented");
}

pair<Mat<ZZ_pX>, long> Operations::pack_rgsw(const Vec<ZZ_pX>& msgs, const long& mode) {
    Mat<ZZ_pX> pack; pack.SetDims(2*l, 2); 

    // clear pack
    for (int i = 0; i < 2*l; i++){
        clear(pack[i][0]);
        clear(pack[i][1]);
    }

    if (mode == 1 || mode == 2) {
        long fator = N / powers[mode-1];
        Mat<ZZ_pX> sum; ZZ_pX power;

        for (long i = 0; i < r; i++) {
            clear(power); SetCoeff(power, i*fator, 1);
            sum = rgsw.enc(msgs[i]);
            for (int j = 0; j< 2*l;j++){
                pack[j][0] += (sum[j][0] * power) % mod_f;
                pack[j][1] += (sum[j][1] * power) % mod_f;
            }
        }
        return {pack, mode};
    }

    if (mode == 3 || mode == 4) {
        long fator1 = N / powers[1 - mode % 2];
        long fator2 = N / powers[mode % 2];
        Mat<ZZ_pX> sum; ZZ_pX power; ZZ_pX db_fator;

        for (long i = 0; i < r; ++i) {
            clear(power); SetCoeff(power, i*fator2, 1);
            db_fator = aut(vec_dbasis[mode - 3][i], fator1, N, mod_f);
            sum = rgsw.enc(msgs[i]);
            for (int j = 0; j< 2*l;j++){
                pack[j][0] += (sum[j][0] * power * db_fator) % mod_f;
                pack[j][1] += (sum[j][1] * power * db_fator) % mod_f;
            }
        }
        return {pack, mode};
    }

    Error("Operation mode not implemented");
}

Vec<ZZ_pX> Operations::unpack_gsw(const pair<Mat<ZZ_pX>, long>& p_gsw) {
    Error("Homomorphic trace in GSW not implemented");
}

ZZ Operations::get_noise(Vec<ZZ_pX> cipher, ZZ_pX msg) {
    return rlwe.get_noise(cipher, msg);  
}   


ZZ Operations::get_noise(Mat<ZZ_pX> cipher, ZZ_pX msg) {
    return rgsw.get_noise(cipher, msg);  
}


ZZ_pX Operations::pack_msgs(Vec<ZZ_pX> msgs, long mode){
    ZZ_pX pack; clear(pack);
    
    if (mode == 1 || mode == 2) {
        long fator = N / powers[mode-1];
        ZZ_pX power;

        for (long i = 0; i < r; i++) {
            clear(power); SetCoeff(power, i*fator, 1);
            pack = pack + (msgs[i]* power) % mod_f;
        }
        return pack;
    }

    if (mode == 3 || mode == 4) {
        long fator1 = N / powers[1 - mode % 2];
        long fator2 = N / powers[mode % 2];

        Vec<ZZ_pX> sum; ZZ_pX power; ZZ_pX db_fator;
        for (long i = 0; i < r; ++i) {
            clear(power); SetCoeff(power, i*fator2, 1);
            db_fator = aut(vec_dbasis[mode-3][i], fator1, N, mod_f);
            pack = pack + (((msgs[i] * power) % mod_f) * db_fator) % mod_f;
        }
        return pack;
    }

    Error("Operation mode not used");
}
