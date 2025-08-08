load("./trace/trace_lib.sage")
load("./GSW/GSW.sage")
load("./CRLWE/RLWE.sage")
load("./key_switch/key_switch_lib.sage")
load("./utils/utils.sage")
load("./utils/gadget_matrix.sage")

Zx = ZZ['x']

def homo_trace_rlwe(ciphertext, P, evk, evk_ktensor, B, q, l, n, prime, exp, base, auto):
    f = Zx(cyclotomic_polynomial(n))
    Rq = (ZZ.quotient(q))['x'].quotient(f)

    a, b = ciphertext[0], ciphertext[1]
    rlwe_c = [Zx(0), a*P.inverse_mod(q)]
    d = ext_prod(B, q, l, n, Rq, rlwe_c, evk)

    k, res = 0, vector(Rq, [0]*2) 
    for i in auto:
        ks_res = ks_rlwe([Rq(Zx(d[0].lift())(x**i)),Rq(Zx(d[1].lift())(x**i))], evk_ktensor[k], base, B, q, n, Rq)
        res = res + ks_res
        k += 1
    
    return vector(Rq, [0, generic_trace(Zx(b.lift()), n, prime, exp)*P.inverse_mod(q)]) - res


def fast_homo_trace_rlwe(rlwe_ciphertext, P, evk, evk_ktensor, B, q, l, n, prime, exp, base, auto):
    f = Zx(cyclotomic_polynomial(n))
    Rq = (ZZ.quotient(q))['x'].quotient(f)

    a, b = rlwe_ciphertext[0], rlwe_ciphertext[1]
    rlwe_c = [Rq(0), a * P.inverse_mod(q)]
    d = ext_prod(B, q, l, n, Rq, rlwe_c, evk)
    
    for j in range(exp):
        res = vector(Rq, [0]*2) 
        k = 0
        for i in auto[j]: 
            ks_res = ks_rlwe([Rq(Zx(d[0].lift())(x**i)),Rq(Zx(d[1].lift())(x**i))], evk_ktensor[j][k], base, B, q, n, Rq)
            res = res + ks_res
            k += 1
        
        d = res     
    return vector(Rq, [0, generic_trace(Zx(b.lift()), n, prime, exp) * P.inverse_mod(q)]) - d

def big_fast_homo_trace_rlwe(rlwe_ciphertext, P, evk, evk_ktensor, B, q, l, n, prime, exp, base, auto):
    f = Zx(cyclotomic_polynomial(n))
    Rq = (ZZ.quotient(q))['x'].quotient(f)

    a, b = rlwe_ciphertext[0], rlwe_ciphertext[1]
    rlwe_c = [Rq(0), a * P.inverse_mod(q)]
    d = ext_prod(B, q, l, n, Rq, rlwe_c, evk)

    for j in range(exp):
        res = vector(Rq, [0]*2) 
        k = 0
        for i in auto[j]:
            ks_res = ks_rlwe([big_aut(d[0], i, Rq),big_aut(d[1], i, Rq)], evk_ktensor[j][k], base, B, q, n, Rq)
            res = res + ks_res
            k += 1
        
        d = res
    
    return vector(Rq, [0, generic_trace(Zx(b.lift()), n, prime, exp) * P.inverse_mod(q)]) - d
