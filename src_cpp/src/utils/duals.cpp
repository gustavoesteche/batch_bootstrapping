#include <NTL/ZZ_pX.h>
#include <NTL/ZZ.h>
#include <iostream>

#include "utils.h"

using namespace std;
using namespace NTL;

Vec<ZZ_pX> canon_dbasis(long p, long n, const ZZ_pXModulus& mod_phi){
    long pn1 = pow(p, n-1);
    long pn = pn1*p;
    long phi_pn = euler_phi(pn);

    Vec<ZZ_pX> dbasis; dbasis.SetLength(phi_pn);
    ZZ_pX power1, power2;
    for (int i=0;i<pn1;i++) {
        clear(dbasis[i]); clear(power1); clear(power2);
        SetCoeff(power1, pn-i, 1);
        power1 = power1 % mod_phi;
        SetCoeff(power2, pn1-i, 1);
        power2 = power2 % mod_phi;
        sub(dbasis[i], power1, power2);
    }

    for (int i=pn1;i<phi_pn;i++) {
        clear(dbasis[i]); clear(power1); clear(power2);
        SetCoeff(power1, pn-i, 1);
        power1 = power1 % mod_phi;
        SetCoeff(power2, pn + pn1-i, 1);
        power2 = power2 % mod_phi;
        sub(dbasis[i], power1, power2); 
        add(dbasis[i], dbasis[i], dbasis[i - pn1]);
    }

    return dbasis;
}