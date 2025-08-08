load("./utils/utils.sage")

Zx = ZZ['x']

def inv_g_ZZ(a, B, q):
    a = sym_mod(ZZ(a), q)
    l = ceil(log(q,B))
    return vector(a.digits(base=B, padto=l))

def inv_g_poly(a, B, q, n):
    l = ceil(log(q,B))
    result = vector(Zx, [0]*l)
    pow_x = Zx(1)
    for i in range(euler_phi(n)):
        result += pow_x * inv_g_ZZ(a[i], B, q)
        pow_x *= Zx(x)
    return result

def inv_g_row_ciphertext(B, q, l, n, Rq, c):
    a, b = c[0], c[1]
    res = vector(Rq, [0] * 2 * l)
    res[0:l] = inv_g_poly(a, B, q, n)
    res[l:2 * l] = inv_g_poly(b, B, q, n)
    
    return res

def ext_prod(B, q, l, n, Rq, rlwe_c, gsw_c):
    return inv_g_row_ciphertext(B, q, l, n, Rq, rlwe_c) * gsw_c

def gadget_vector(B, q):
    l = ceil(log(q,B))
    return vector([B**i for i in range(l)])