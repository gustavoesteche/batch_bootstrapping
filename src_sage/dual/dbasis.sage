load("./trace/trace_lib.sage")

Z = PolynomialRing(ZZ, 'x')
R = PolynomialRing(QQ, 'x') 

# In the toolkit article the dual can be defined as < 1 - zeta_p / p^m > = < 1 - x^(p^(m-1)) / p^m >
def dual_factor(p, m):
    '''Returns the dual generator, in which the elements of the dual must be ring multiples'''
    return Qx((1 - x ** (p ** (m-1))) / (p**m)) % Qx(cyclotomic_polynomial(p ** m))
    
def sample_dual(p, m):
    '''Sample a random element from the dual'''
    return (Qx(Zx.random_element()) * dual_factor(p,m)) % Qx(cyclotomic_polynomial(p ** m))

# The strategy of computing t
def dual_trace_pm_to_pm_1(dual_element, pd, md):
    '''Compute the trace from the field pm to pm-1'''
    return fast_trace_pm_to_pm_1(pd, md, dual_element * (pd ** md)) / (pd ** md)

def dual_trace_pm_to_1(dual_element, pd, md):
    '''Compute the trace from the field pm to 1'''
    return fast_trace_pm_to_1(pd, md, dual_element * (pd ** md)) / (pd ** md) 

def deb_canon_dbasis(p, m, dbasis, basis):
    phi_pm = euler_phi(p**m)
    fm = R(cyclotomic_polynomial(p**m)) 
    
    for i in range(phi_pm):
        print(i, dbasis[i])
        for j in range(phi_pm):
            print(trivial_trace_m_to_1(dbasis[i] * basis[j]*p**m, p**m))

Zx = PolynomialRing(ZZ, 'x')
Qx = PolynomialRing(QQ, 'x') 

def print_traces(p, m):
    pm = p **m
    phi_pm = euler_phi(p ** m)
    for i in range(phi_pm):
        print(fast_trace_pm_to_1(p, m, Zx(x^i)))

def check_dual_basis(basis, dual_basis, p, n):
    m = p**n
    for i in range(len(basis)):
        for j in range(len(dual_basis)):
            if i != j:
                assert fast_trace_pm_to_1(p, n, basis[i] * (dual_basis[j] * m)) == 0, "basis[{}]:{} dbasis[{}]:{}".format(i, basis[i], j, dual_basis[j])
            else: 
                assert fast_trace_pm_to_1(p, n , basis[i] * (dual_basis[j]*m))/m == 1,"basis[{}]:{} dbasis[{}]:{}".format(i, basis[i], j, dual_basis[j])

def canon_dbasis(p, m):
    f = Zx(cyclotomic_polynomial(p**m))
    l = vector(Qx, [0]*euler_phi(p**m))
    for i in range(p**(m-1)):
        p_magic = (p**(m-1) - i) % (p**m)
        fj = (Zx(x**(p**m - i) - x**p_magic) % f)
        l[i] = fj
    
    for i in range(p**(m-1), euler_phi(p**m)):
        p_magic = (p**(m-1) - i) % (p**m)
        fj = (Zx(x**(p**m - i) - x**p_magic) % f)
        fj = fj + l[i - p**(m-1)]
        l[i] = fj
        
    return l/(p**m)
