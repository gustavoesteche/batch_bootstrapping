#ifndef UTILS_H
#define UTILS_H

#include <NTL/ZZ_pEX.h>
#include <NTL/ZZX.h>
#include <NTL/vector.h>
#include <iostream>

using namespace std;
using namespace NTL;

/**
 * @brief Converts an integer-coefficient polynomial to a modular polynomial.
 * @param f Input polynomial over ZZ.
 * @return Polynomial over ZZ_pX with coefficients reduced modulo current modulus.
 */
ZZ_pX conversion_X_to_pX(const ZZX &f);

/**
 * @brief Converts a modular polynomial to an integer-coefficient polynomial.
 * @param f Input polynomial over ZZ_pX.
 * @return Polynomial over ZZ with lifted coefficients.
 */
ZZX conversion_pX_to_X(const ZZ_pX &f);

/**
 * @brief Computes Euler's totient function φ(n).
 * @param n Positive integer.
 * @return φ(n).
 */
long euler_phi(long n);

/**
 * @brief Generates the N-th cyclotomic polynomial over ZZ.
 * @param N Cyclotomic order.
 * @return N-th cyclotomic polynomial.
 */
ZZX cyclotomic_polynomial(long N);

/**
 * @brief Multiplies a vector of polynomials by a polynomial matrix modulo mod_phi.
 * @param res Output vector.
 * @param v Input polynomial vector.
 * @param M Input polynomial matrix.
 * @param mod_phi Precomputed modulus for modular polynomial arithmetic.
 */
void Mul(Vec<ZZ_pX> &res, const Vec<ZZ_pX>& v, const Mat<ZZ_pX>& M, const ZZ_pXModulus& mod_phi);

/**
 * @brief Applies an automorphism to a polynomial: x ↦ x^alpha mod Φ_N(x).
 * @param poly Input polynomial.
 * @param alpha Automorphism exponent.
 * @param N Cyclotomic order.
 * @param mod_phi Precomputed modulus for Φ_N(x).
 * @return Transformed polynomial.
 */
ZZ_pX aut(const ZZ_pX& poly, long alpha, const long N, const ZZ_pXModulus& mod_phi);

/**
 * @brief Finds automorphism exponents for a given cyclotomic order and prime power.
 * @param N Cyclotomic order.
 * @param p Prime base.
 * @param exp Exponent.
 * @return Vector of automorphism exponents.
 */
Vec<long> find_aut(long N, long p, long exp);

/**
 * @brief Finds automorphism exponents for a given cyclotomic order and prime.
 * @param N Cyclotomic order.
 * @param p Prime base.
 * @return Vector of automorphism exponents.
 */
Vec<long> find_aut_tower(long N, long p);

/**
 * @brief Constructs powers of automorphisms for fast evaluation.
 * @param N Cyclotomic order.
 * @param p Prime base.
 * @param exp Exponent.
 * @return Table of automorphism powers.
 */
Vec<Vec<long>> construct_aut_power(long N, long p, long exp);

/**
 * @brief Reduces each coefficient of a polynomial modulo a.
 * @param poly Input polynomial over ZZ_pX.
 * @param a Modulus for coefficient reduction.
 * @return Reduced polynomial.
 */
ZZ_pX mod_poly(const ZZ_pX &poly, long a);

/**
 * @brief Reduces each coefficient of a polynomial modulo a.
 * @param poly Input polynomial over ZZ.
 * @param a Modulus for coefficient reduction.
 * @return Reduced polynomial.
 */
ZZX mod_poly(const ZZX &poly, long a);

/**
 * @brief Applies centered modular reduction to an integer.
 * @param a Integer to reduce (in-place).
 * @param q Modulus.
 */
void sym_mod(ZZ &a, ZZ q);

/**
 * @brief Applies centered modular reduction to a polynomial.
 * @param poly Polynomial over ZZ_pX.
 * @param q Modulus.
 * @return Polynomial with coefficients in (-q/2, q/2].
 */
ZZX sym_mod_poly(const ZZ_pX &poly, ZZ q);

/**
 * @brief Applies centered modular reduction to a vector of polynomials.
 * @param vec Vector of polynomials over ZZ_pX.
 * @param q Modulus.
 * @return Vector of reduced polynomials.
 */
Vec<ZZX> sym_mod_vec(const Vec<ZZ_pX>& vec, ZZ q);

/**
 * @brief Computes the infinity norm (max absolute coefficient) of a polynomial.
 * @param poly Polynomial over ZZ.
 * @return Infinity norm as ZZ.
 */
ZZ infinity_norm(ZZX poly);

/**
 * @brief Computes the infinity norm of a vector of polynomials.
 * @param vec Vector of polynomials over ZZ.
 * @return Infinity norm as ZZ.
 */
ZZ infinity_norm_vec(Vec<ZZX> vec);

/**
 * @brief Inverse gadget decomposition of a scalar.
 * @param a_in Input modular integer.
 * @param q Modulus.
 * @param l Gadget decomposition length.
 * @return Vector of decomposed coefficients over ZZ.
 */
Vec<ZZ> inv_g_ZZ(const ZZ_p &a_in, const ZZ& q, long l);

/**
 * @brief Inverse gadget decomposition of a polynomial over ZZ_pX to ZZ.
 * @param a Input polynomial.
 * @param q Modulus.
 * @param l Gadget decomposition length.
 * @param n Polynomial degree.
 * @return Vector of decomposed polynomials over ZZ.
 */
Vec<ZZX> inv_g_poly(const ZZ_pX& a, const ZZ& q, long l, long n);

/**
 * @brief Inverse gadget decomposition of a polynomial over ZZ_pX.
 * @param a Input polynomial.
 * @param q Modulus.
 * @param l Gadget decomposition length.
 * @param n Polynomial degree.
 * @return Vector of decomposed polynomials over ZZ_pX.
 */
Vec<ZZ_pX> inv_g_poly_p(const ZZ_pX& a, const ZZ& q, long l, long n);

/**
 * @brief Inverse gadget decomposition of a RLWE ciphertext.
 * @param a Input ciphertext polynomials.
 * @param q Modulus.
 * @param l Gadget decomposition length.
 * @param n Polynomial degree.
 * @return Vector of decomposed ciphertexts over ZZ_pX.
 */
Vec<ZZ_pX> inv_g_ciphertext_p(const Vec<ZZ_pX>& a, const ZZ& q, long l, long n);

/**
 * @brief Computes the external product between a RLWE ciphertext and a RGSW ciphertext.
 * @param rlwe_c RLWE ciphertext (vector of polynomials).
 * @param gsw_c RGSW ciphertext (matrix of polynomials).
 * @param q Modulus.
 * @param l Gadget decomposition length.
 * @param n Polynomial degree.
 * @param mod_phi Precomputed modulus for Φ_N(x).
 * @return RLWE ciphertext after external product.
 */
Vec<ZZ_pX> ext_prod(const Vec<ZZ_pX>& rlwe_c, Mat<ZZ_pX> gsw_c, const ZZ& q, long l, long n, const ZZ_pXModulus &mod_phi);

#endif