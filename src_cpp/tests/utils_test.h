#ifndef UTILS_TEST_H
#define UTILS_TEST_H

#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <iostream>

using namespace NTL;
using namespace std;

#include "../src/utils/utils.h"

void build_poly(ZZX &f, ZZ_pXModulus &mod_f, const long N); 

Vec<pair<long,long>> build_prime_power(pair<long, long> p1, pair<long, long> p2, pair<long, long> p3);

ZZ_pX generate_message(long N, const ZZ_pXModulus& mod_phi);

ZZ_pX generate_message_depth(long N);

Vec<ZZ_pX> sample_packages(long N, long r, Vec<pair<long, long>> pp, const ZZ_pXModulus& mod_phi, bool depth=false);

bool compare_poly(const ZZ_pX& p1, const ZZ_pX& p2);

bool compare_poly(const ZZX& p1, const ZZX& p2);

bool compare_vec(const Vec<ZZ_pX>& v1, const Vec<ZZ_pX>& v2);

#endif