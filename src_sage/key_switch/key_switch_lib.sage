load("./utils/gadget_matrix.sage")
load("./utils/utils.sage")

Zx = PolynomialRing(ZZ, "x")

def ks_rlwe(ciphertext, K, base, B, q, N, Rq):
    a, b = ciphertext[0], ciphertext[1]
    ev_aK = inv_g_poly(a, base, q, N) * K
    ev_aK = [Rq(ev_aK[0]), Rq(ev_aK[1])]
    
    return vector([-ev_aK[0], b - ev_aK[1]])
    
def ks_gsw(ciphertext, ks_key, K, base, B, q, N, Rq): 
    l = ceil(log(q, B))
    new_C = Matrix(Rq, 2* l, 2)
    
    for i in range(l):
        new_C[l+i] = ks_rlwe(ciphertext[l+i], K, base, B, q, N, Rq)
        new_C[i] = inv_g_row_ciphertext(B, q, l, n, Rq, new_C[l+i]) * ks_key
    return new_C
