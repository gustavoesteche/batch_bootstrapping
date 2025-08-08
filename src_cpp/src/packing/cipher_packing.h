#ifndef CIPHER_PACKING_H
#define CIPHER_PACKING_H 

#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_p.h>
#include <iostream>
#include <vector>

#include "../utils/utils.h"
#include "../utils/trace.h"
#include "../utils/duals.h"
#include "../utils/key_switch.h"
#include "../utils/RLWE_distribution.h"
#include "../schemes/RLWE.h"
#include "../schemes/RGSW.h"

using namespace std;
using namespace NTL;

class Operations {
public:
    Operations(const long &N, const ZZ& modulus, Vec<pair<long, long>> pp, double sigma = 3.2);

    /**
     * @brief Computes the trace of a RLWE ciphertext for a given mode/index.
     * @param poly RLWE ciphertext
     * @param index Target index for the trace.
     * @return Resulting RLWE ciphertext vector after trace.
     */
    Vec<ZZ_pX> trace(Vec<ZZ_pX> poly, long index);

    /**
     * @brief Performs the external product described in this framework.
     * @param p_rlwe RLWE ciphertext (vector of polynomials, mode).
     * @param p_gsw  RGSW ciphertext (matrix of polynomials, mode).
     * @return RLWE ciphertext after external product and trace.
     */
    pair<Vec<ZZ_pX>, long> ext_prod_trace(pair<Vec<ZZ_pX>, long> p_rlwe, pair<Mat<ZZ_pX>, long> p_gsw);

    /**
     * @brief Packs a vector of messages into a single RLWE ciphertext.
     * @param msgs Messages to pack.
     * @param mode Packing mode.
     * @return Packed RLWE ciphertext (vector of polynomials, mode).
     */
    pair<Vec<ZZ_pX>, long> pack_rlwe(const Vec<ZZ_pX> & msgs, const long& mode);

    /**
     * @brief Unpacks a RLWE ciphertext into individual messages.
     * @param p_rlwe Packed RLWE ciphertext.
     * @return Vector of unpacked message polynomials.
     */
    Vec<ZZ_pX> unpack_rlwe(const pair<Vec<ZZ_pX>, long>& p_rlwe);

    /**
     * @brief Packs a vector of messages into a RGSW ciphertext.
     * @param msgs Messages to pack.
     * @param mode Packing mode.
     * @return Packed RGSW ciphertext (matrix of polynomials, mode).
     */
    pair<Mat<ZZ_pX>, long> pack_rgsw(const Vec<ZZ_pX>& msgs, const long& mode);

    /**
     * @brief Unpacks a RGSW ciphertext into individual messages.
     * @param p_gsw Packed RGSW ciphertext.
     * @return Vector of unpacked message polynomials.
     */
    Vec<ZZ_pX> unpack_gsw(const pair<Mat<ZZ_pX>, long>& p_gsw);

    /**
     * @brief Computes the infinite norm of the noise noise in a RLWE ciphertext.
     * @param cipher RLWE ciphertext.
     * @param msg Expected plaintext polynomial.
     * @return Infinity norm of the noise as a ZZ integer.
     */
    ZZ get_noise(Vec<ZZ_pX> cipher, ZZ_pX msg);

    /**
     * @brief Computes the infinite norm of the noise in a RGSW ciphertext.
     * @param cipher RGSW ciphertext.
     * @param msg Expected plaintext polynomial.
     * @return Infinity norm of the noise as a ZZ integer.
     */
    ZZ get_noise(Mat<ZZ_pX> cipher, ZZ_pX msg);
    
    /**
     * @brief Packs multiple message polynomials into a single polynomial.
     * @param msgs Messages to pack.
     * @param mode Packing mode.
     * @return Packed polynomial.
     */
    ZZ_pX pack_msgs(Vec<ZZ_pX> msgs, long mode);

    
private:    
    long N, l, r, phi_N;
    ZZ Q;
    Vec<pair<long, long>> pp;
    vector<long> powers;  
    
    RLWE rlwe; 
    RGSW rgsw;
    
    vector<Vec<Vec<long>>> auto_ks;
    vector<Vec<Mat<ZZ_pX>>> evk_ktensor;
    vector<Vec<ZZ_pX>> vec_dbasis;
    Mat<ZZ_pX> evk;
      
    ZZX f;
    ZZ_pXModulus mod_f;

    void construct_aut();
    void gen_eval_keys();
    void gen_canon_dbasis();
};

#endif
