load("./tensor_composition/trace_composition.sage")

Zx = PolynomialRing(ZZ, 'x')

# testing for the composition of twwo rings
def trace_r1r2_to_r2(poly1, p1, m1, poly2):
    return Zx(fast_trace_pm_to_1(p1, m1, poly1) * poly2)

def trace_Rr1_to_R(poly1, p1, m1, poly_m, index = 0):
    '''
    Computes Tr_{Rr1/R} 
    index: which polynomial will be multiplied by the trace factor
    poly1: the polynomial in the field r1 thats going to be left out=
    p1, m1: the prime powers that defines the ring of integers r1
    poly_m: the decomposed polynomials of R
    '''
    trace = fast_trace_pm_to_1(p1, m1, poly1)

    trace_res = []
    for poly in poly_m:
        trace_res.append(poly)

    trace_res[index] = trace_res[index] * trace 
    return trace_res

# Calculating Trace R1 x R2 x R3 / R1 x R(2|3) 
def trace_mode(poly_m, factors_m, mode, index = 0):
    '''
    Computes Tr_K/K13 mode = 1 or Tr_K/K12 mode = 2
    index: which polynomial will be multiplied by the trace factor
    factors_m: A tuple list representing the prime power that divides m
    poly_m[i]: the polynomial in the field K_{p_i ^m_i} that decomposes the original poly 
    '''
    assert(len(poly_m) == 3), "3 fields are necessary"
    assert(mode == 1 or mode == 2), "Choose mode 1 or 2"

    poly_n = []
    for i in range(len(poly_m)):
        if i != mode:  
            poly_n.append(poly_m[i])
    
    return trace_Rr1_to_R(poly_m[mode], factors_m[mode][0],factors_m[mode][1], poly_n, index)


# Test Trace Tr_{R1xR2/R1}
def test_1():
    p0, m0 = next_prime(ZZ.random_element(1, 11)) , ZZ.random_element(1, 3)
    while true:
        p1, m1 = next_prime(ZZ.random_element(1, 11)) , ZZ.random_element(1, 3)
        if (p1 != p0):
            break
    
    poly0 = Zx.random_element(degree= euler_phi(p0 ** m0)-1)

    polym = [poly0, poly1]
    fm = [(p0,m0),(p1,m1)]

    comp_trace = composed_trace(fm, polym)
    trace1 = trace_r1r2_to_r2(poly0, p0, m0, poly1)
    trace2 = fast_trace_pm_to_1(p1, m1, trace1)

    trace_1 = trace_r1r2_to_r2(poly1, p1, m1, poly0)
    trace_2 = fast_trace_pm_to_1(p0, m0, trace_1)

    #print(comp_trace)
    #print(trace2)
    #print(trace_2)
    assert (comp_trace == trace2 == trace_2)

# Test Trace Tr_{R1xR2xR3/R1xR2}
def test_2():
    p0, m0 = 3 , 2
    p1, m1 = 5 , 2
    p2, m2 = 7 , 2
    
    poly0 = Zx.random_element(degree = euler_phi(p0 ** m0)-1)
    poly1 = Zx.random_element(degree = euler_phi(p1 ** m1)-1)
    poly2 = Zx.random_element(degree = euler_phi(p2 ** m2)-1)
    poly_m01 = [poly0, poly1]
    poly_m = [poly0, poly1, poly2]

    fm01 = [(p0,m0),(p1,m1)]
    fm = [(p0,m0),(p1,m1),(p2,m2)]

    comp_trace = composed_trace(fm, poly_m)

    trace1 = trace_Rr1_to_R(poly2, p2, m2, poly_m01)
    trace2 = composed_trace(fm01, trace1)

    trace_1 = trace_mode(poly_m, fm, 2)
    trace_2 = composed_trace(fm01, trace_1)

    #print(comp_trace)
    #print(trace2)
    #print(trace_2)
    assert (comp_trace == trace2 == trace_2)

for i in range(10):
    #print(i)
    test_2()
    
print("traces between the tensor are working")