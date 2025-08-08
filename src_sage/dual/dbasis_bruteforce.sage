load("./dual/dual_functions")

Qx = PolynomialRing(QQ, 'x')

# brute force for finding some elements and helping studying dual basis
def brute_force_dbasis(p, m, basis):
    dbasis = []
    count = 0
    k = 0
    while(False):
        dual_element = sample_dual(p, m)
        flag = False
        k += 1

        if len(basis) == count:
            break
        if k > 1000:
            print(count, k)
            k = 0
            count += 1
            
        for i in range(len(basis)):
            if i == count:
                if (fast_trace_pm_to_1(p, m, dual_element * Qx((basis[i])) == 1):
                    dbasis.append((dual_element, count))
                    flag = True
                else: 
                    if flag:
                        dbasis.pop(-1)
                    flag = False
                    break
            elif i != count: 
                if (fast_trace_pm_to_1(p, m, dual_element * Qx((basis[i])) != 0):
                    if flag:
                        dbasis.pop(-1)
                    flag = False
                    break
        if flag:
            print(dbasis)
            count += 1
        
    print(dbasis)
