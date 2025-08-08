#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz.h>

#include "cyclotomic.h"

void cyclotomic_coeffs(long N, long* length, long *result) {
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpz_poly_cyclotomic(poly, N);
    *length = fmpz_poly_length(poly);

    if (!result) {
        fmpz_poly_clear(poly);
        return NULL;
    }

    for (long i = 0; i < *length; i++) {
        result[i] = fmpz_get_si(fmpz_poly_get_coeff_ptr(poly, i));
    }

    fmpz_poly_clear(poly);
}
