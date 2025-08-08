#ifndef RLWE_DISTRIBUTION_H
#define RLWE_DISTRIBUTION_H

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <random>
#include <cmath>

using namespace std;
using namespace NTL;

class DiscreteGaussian {
public:
    DiscreteGaussian() = default;
    DiscreteGaussian(double sigma);
    ZZ sample();
private:
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
};

class RLWE_Distribution {
public:
    RLWE_Distribution() = default;
    RLWE_Distribution(const ZZ_pX& sk, long N, const ZZ& modulus, double sigma = 3.2);

    NTL::ZZ_pX random_noise();
    NTL::ZZ_pX random_a();
    NTL::Vec<NTL::ZZ_pX> sample();

private:
    long n;
    long phi_n;
    ZZ q;
    ZZ_pX s;
    ZZX f;
    DiscreteGaussian D;
    ZZ_pXModulus mod_f;
};

#endif
