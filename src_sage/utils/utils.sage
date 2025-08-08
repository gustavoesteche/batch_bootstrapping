def sym_mod(a, q):
    a = ZZ(a) % q
    if 2*a > q:
        return a - q
    return ZZ(a)

def sym_mod_poly(poly, q):
    return Zx([sym_mod(ZZ(ai), q) for ai in poly.list()])

def sym_mod_vec(vec, q):
    return [sym_mod_poly(vi, q) for vi in vec]

def round_poly(h):  
    return Zx([round(hi) for hi in h.list()])

def infinity_norm(poly):
    if 0 == poly:
        return 0
    return vector(ZZ, poly.list()).norm(Infinity)

def infinity_norm_vec(vec):
    return max([infinity_norm(vi) for vi in vec])

def aut(poly, m, i):
    return Zx(poly(x=x^i)) % Zx(cyclotomic_polynomial(m))

def big_aut(poly, alfa, R):
    res = R(0)
    for i in range(len(list(poly))):
        res = res + poly[i] * R(x**(i*alfa)) 
    return res

def find_aut(m, p, n):
    mi = p ** n
    fator = m/mi
    
    k = p-(fator % p).inverse_mod(p)
    rest = [i for i in range(p) if i!=k]
    p_power = mi//p
    aut = [(i*p + j)*fator+1 for i in range(p_power) for j in rest]
    return aut

def find_aut_tower(m, p):
    fator = m/p
    if fator % p == 0:
        return [(i*fator+1) for i in range(p)]
    else: 
        return find_aut(m, p, 1)

def construct_aut_tower(m, p):
    auto = []
    while(m % p == 0):
        auto.append(find_aut_tower(m, p))
        m = m // p 
    return auto

