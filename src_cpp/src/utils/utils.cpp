#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <iostream>

extern "C" {
    #include "cyclotomic.h"
}

using namespace std;
using namespace NTL;

ZZ_pX conversion_X_to_pX(const ZZX &f) {
    ZZ_pX g;
    long deg_ = deg(f);
    g.SetMaxLength(deg_ + 1);

    for (long i = 0; i <= deg_; ++i) {
        ZZ coeff_ = coeff(f, i);
        SetCoeff(g, i, to_ZZ_p(coeff_));
    }

    return g;
}

ZZX conversion_pX_to_X(const ZZ_pX &f) {
    ZZX g;
    long deg_ = deg(f);

    for (long i = 0; i <= deg_; ++i) {
        SetCoeff(g, i, rep(coeff(f, i)));  // rep extracts the ZZ from ZZ_p
    }

    return g;
}

ZZ_pX mod_poly(const ZZ_pX &poly, long a) {
    ZZ_pX ans; clear(ans);
    for (int i = 0; i<= deg(poly); i++) {
        SetCoeff(ans, i, rem(ZZ(conv<long>(poly[i])), a));
    } 
    return ans;
}

ZZX mod_poly(const ZZX &poly, long a) {
    ZZX ans; clear(ans);
    for (int i = 0; i<= deg(poly); i++) {
        SetCoeff(ans, i, rem(ZZ(conv<long>(poly[i])), a));
    } 
    return ans;
}

void Mul(Vec<ZZ_pX> &res, const Vec<ZZ_pX>& v, const Mat<ZZ_pX>& M, const ZZ_pXModulus& mod_phi) {
    long m = M.NumCols();
    long n = M.NumRows();

    res.SetLength(m);
    for (long j = 0; j < m; ++j) {
        clear(res[j]);
        for (long i = 0; i < n; ++i){    
            res[j] += v[i] * M[i][j];
            rem(res[j], res[j], mod_phi);
        }
    }
}

long euler_phi(long n) {
    long result = n; 
    
    for (long long p = 2; p * p <= n; ++p) {
        if (n % p == 0) {
            result -= result / p;
            while (n % p == 0) {
                n /= p;
            }
        }
    }

    if (n > 1) {
        result -= result / n;
    }
    
    return result;
}

ZZX cyclotomic_polynomial(const long N) {
    long result[N];
    long* length = new long;
    cyclotomic_coeffs(N, length, result);
    ZZX f; clear(f);
    for (int i=0;i<*length;i++) {
        SetCoeff(f, i, result[i]);
    }

    delete length;
    return f;
}


ZZ_pX aut(const ZZ_pX& poly, long alpha,  const long N, const ZZ_pXModulus& mod_phi) {    
    ZZ_pX res; clear(res);

    ZZ_pX term, power;
    for (long i = 0; i <= deg(poly); ++i) {
        if (IsZero(poly[i])) continue;
        SetX(power); 
        PowerMod(power, power, i * alpha % N, mod_phi);
        mul(term, power, poly[i]); 
        add(res, res, term);
    }

    return res;
}

ZZ_pX apply_automorphism(const ZZ_pX& f, long alpha, const ZZ_pXModulus& mod_phi) {
    ZZ_pX result;
    
    ZZ_pX xk;
    power(xk, ZZ_pX::zero()+1, alpha); 

    CompMod(result, f, xk, mod_phi);
    return result;
}

Vec<long> find_aut(long N, long prime, long exp){
    long mi1; mi1 = pow(prime, exp-1);
    long mi = mi1 * prime;
    long fator = N / mi; 
    long k = prime - conv<long> ( InvMod(fator, prime) );

    long q = 0;
    Vec<long> res; res.SetLength(mi1*(prime-1)); 
    for (int i=0;i<mi1;i++) {
        for (int j=0;j<prime;j++) {
            if(j!=k)
                res[q++] = (i*prime+j)*conv<long>(fator)+1;
        }
    }
    return res;
}

Vec<long> find_aut_tower(long N, long prime) {
    long fator = N / prime;
    if (fator % prime == 0) {
        Vec<long> res; res.SetLength(prime);
        for (int i = 0; i < prime; i++) 
            res[i] = i*fator+1;
        return res;
    }
    return find_aut(N,prime,1);
}

Vec < Vec<long> > construct_aut_power(long N, long p, long exp) {
    Vec<Vec<long>> res; res.SetLength(exp);
    Vec<long> loc;
    for (int i = 0; i<exp; i++) {
        loc = find_aut_tower(N, p);
        res[i].SetLength(loc.length());
        for (int j = 0; j<loc.length(); j++) {
            res[i][j] = loc[j];
        }
        N /= p;
    }
    return res;
}

void sym_mod(ZZ &a, ZZ q){
    a %= q;
    if(2* a > q)
        a -= q;
}

ZZ sym_mod(const ZZ_p &a, ZZ q){
    ZZ new_a = conv<ZZ>(a) % q;
    if(2* new_a > q){
        new_a -= q;
    }
    return new_a;
}

ZZX sym_mod_poly(const ZZ_pX &poly, ZZ q){
    ZZ t; 
    ZZX new_poly; clear(new_poly);

    for (long i = 0; i <= deg(poly); ++i) {
        SetCoeff(new_poly, i, sym_mod(poly[i], q));
    }

    return new_poly;
}

Vec<ZZX> sym_mod_vec(const Vec<ZZ_pX>& vec, ZZ q){
    Vec<ZZX> new_vec; new_vec.SetLength(vec.length());

    for (long i = 0; i < vec.length(); i++) {
        new_vec[i] = sym_mod_poly(vec[i], q);
    }

    return new_vec;
}  

ZZ infinity_norm(ZZX poly) {
    ZZ max = ZZ(0);
    for (long i = 0; i <= deg(poly); ++i) {
        ZZ a = abs(coeff(poly, i));
        if (a > max) max = a;
    }
    return max;
}

ZZ infinity_norm_vec(Vec<ZZX> vec){
    ZZ maxi = ZZ();
    for (int i=0;i<vec.length();i++) {
        maxi = (maxi < infinity_norm(vec[i])) ? infinity_norm(vec[i]) : maxi;
    }
    return maxi;
}

Vec<ZZ> inv_g_ZZ(const ZZ_p &a_in, const ZZ& q, long l) {
    ZZ a = sym_mod(a_in, q); ZZ sign = ZZ(1);
    Vec<ZZ> res; res.SetLength(l);
    
    if (a < ZZ(0))
        sign = -1;

    for (long i = 0; i < l; ++i) {
        res[i] = sign * bit(a, i);
    }
    return res;
}

Vec<ZZX> inv_g_poly(const ZZ_pX& a, const ZZ& q, long l, long n) {
    long phi_n = euler_phi(n);

    Vec<ZZX> result;
    result.SetLength(l);
    for (long i = 0; i < l; ++i) {
        result[i].SetMaxLength(phi_n);
        clear(result[i]);
    }

    for (long i = 0; i < phi_n; ++i) {
        Vec<ZZ> digits = inv_g_ZZ(coeff(a, i), q, l);
        for (long j = 0; j < l; ++j) {
            SetCoeff(result[j], i, coeff(result[j], i) + digits[j]);
        }
    }

    return result;
}

Vec<ZZ_pX> inv_g_poly_p(const ZZ_pX& a, const ZZ& q, long l, long n) {
    long phi_n = euler_phi(n); 

    Vec<ZZ_pX> result;
    result.SetLength(l);
    for (long i = 0; i < l; ++i) {
        result[i].SetMaxLength(phi_n);
        clear(result[i]);
    }

    for (long i = 0; i < phi_n; ++i) {
        Vec<ZZ> digits = inv_g_ZZ(coeff(a, i), q, l);
        for (long j = 0; j < l; ++j) {
            SetCoeff(result[j], i, coeff(result[j], i) + to_ZZ_p(digits[j]));
        }
    }

    return result;
}

Vec<ZZ_pX> inv_g_ciphertext_p(const Vec<ZZ_pX>& a, const ZZ& q, long l, long n){
    Vec<ZZ_pX> res; res.SetLength(2*l);
    for (int i = 0; i < 2; i++) {
        Vec<ZZ_pX> ans = inv_g_poly_p(a[i], q, l, n);
        for (int j = 0; j<l; j++) {
            res[i*l+j] = ans[j];
        }
    }
    return res;
}

Vec<ZZ_pX> ext_prod(const Vec<ZZ_pX>& rlwe_c, Mat<ZZ_pX> gsw_c, const ZZ& q, long l, long n, const ZZ_pXModulus &mod_phi){
    Vec<ZZ_pX> res; 
    Mul(res, inv_g_ciphertext_p(rlwe_c, q, l, n), gsw_c, mod_phi);
    return res;
}
