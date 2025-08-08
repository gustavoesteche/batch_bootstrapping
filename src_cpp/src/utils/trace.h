#ifndef TRACE_H
#define TRACE_H

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <iostream>

using namespace std;
using namespace NTL;

ZZ_pX generic_trace(const ZZ_pX& poly, long m, long p, long n, const ZZ_pXModulus& mod_phi);

Vec<ZZ_pX> homo_generic_trace(const Vec<ZZ_pX>& rlwe_cipher, ZZ_p p_inverse, Mat<ZZ_pX> evk, Vec<Mat<ZZ_pX>> K_tower, 
    const ZZ Q, const long l, const long N, const long prime, const long exp, Vec< Vec<long> > autos, const ZZ_pXModulus& mod_phi);

#endif