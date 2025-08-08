load("./trace/trace_lib.sage")
load("./tensor_composition/trace_composition.sage")

Zx = PolynomialRing(ZZ, 'x')
Qx = PolynomialRing(QQ, 'x')
 

def dual_factor(p, m):
    '''Returns the dual generator, in which the elements of the dual must be integer multiples'''
    return Qx((1 - x ** (p ** (m-1))) / (p**m))

def sample_dual(p, m):
    '''Sample a random element from the dual'''
    return (Qx(Zx.random_element()) * dual_factor(p,m)) % Qx(cyclotomic_polynomial(p ** m))

def compose_trace_wdual(poly_m, factors_m, dual_element, pd, md):
    '''Compose the trace with a dual field'''
    p_poly_m = poly_m.copy()
    f_factors_m = factors_m.copy()

    p_poly_m.append(dual_element * (pd ** md))
    f_factors_m.append((pd, md))
    
    return compose_trace(f_factors_m, p_poly_m) / (pd ** md)
    
def dual_trace_pm_to_pm_1(dual_element, pd, md):
    '''Compute the trace from the field pm to pm-1'''
    return fast_trace_pm_to_pm_1(pd, md, dual_element * (pd ** md)) / (pd ** md)

def dual_trace_pm_to_1(dual_element, pd, md):
    '''Compute the trace from the field pm to 1'''
    return fast_trace_pm_to_1(pd, md, dual_element * (pd ** md)) / (pd ** md) 

'''
p1,m1 = 5, 3
p2,m2 = 7, 2
p3,m3 = 2, 3
 
factors_m = [(p1,m1),(p2,m2)]
poly_m = [Zx.random_element(degree = euler_phi(p1**m1)-1), Zx.random_element(degree = euler_phi(p2**m2)-1)]
dual_element = sample_dual(p3, m3)
element = Zx.random_element(degree = euler_phi(p3**m3)-1)

print("dual element ",dual_element)
print(dual_element * Qx(element))

# print("trace dual p^m to p^(m-1)",dual_trace_pm_to_pm_1(dual_element, p3,m3))
print("trace dual p^m to 1",fast_trace_pm_to_1(p3, m3, dual_element * Qx(element)))

print("prime powers ", factors_m, (p3, m3))
print("poly composition ", poly_m)
print("Composed dual trace ", compose_trace_wdual(poly_m, factors_m, dual_element, p3, m3))
'''