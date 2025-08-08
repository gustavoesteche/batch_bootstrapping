#ifndef DUAL_H
#define DUAL_H

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <iostream>

#include "utils.h"

using namespace std;
using namespace NTL;

/**
 * @brief Calculates the dual basis of the power basis $p^n$
 * @param p prime
 * @param n exponent
 * @param mod_phi the module of the cyclotomic ring
 * @return Dual basis vector 
 */
Vec<ZZ_pX> canon_dbasis(long p, long n, const ZZ_pXModulus& mod_phi);

#endif