from time import perf_counter 

load("./trace/trace_lib.sage")

def test_pm_to_1():
    ''' ''' 
    def test(p, m, test_size): 
        assert m < 30, "since the operation is exp(m), this size of m would take too long for testing"
        euler_phi_pm = (p-1) * p^(m-1)

        start_time = perf_counter()
        for test_number in range(test_size):
            # 
            poly =Zx.random_element(euler_phi_pm - 1)

            # Computing the trace using the three known functions  
            poly_trivial = trivial_trace_m_to_1(poly, p ** m)
            poly_fast = fast_trace_pm_to_1(p, m, poly)
            poly_calculated = trace_pm_to_1(p, m, poly)   

            # assert the equality 
            assert poly_calculated == poly_fast == poly_trivial, "Test {}/{} Failed \n {}  \n calculated {} \n fast calculated{}".format(test_number+1, test_size, poly_trivial, poly_calculated, poly_fast)
    
        end_time  = perf_counter()
        #print("Test completed, time elapsed: {:2f}".format(end_time - start_time))
    
    def test_multiple(p_range, m_max, tests = 10, test_size = 10):
        #print("\n######\n Testing trace_pm_to_1 \n######\n")
        for i in range(tests):
            p = random_prime(p_range)
            m = ZZ.random_element(1, m_max)
            #print("Teste {}:".format(i+1))
            #print("Prime p = {}, exponent m = {}".format(p, m))
            test(p, m, test_size)
        print("trace pm to 1 working as intended")
    
    test_multiple(p_range = 5, m_max = 4, tests = 5, test_size = 2)

def test_pm_to_pm_1():
    ''' ''' 
    def test(p, m, test_size): 
        assert m < 30, "since the operation is exp(m), this size of m would take too long for testing"
        euler_phi_pm = (p-1) * p^(m-1)

        start_time = perf_counter()
        for test_number in range(test_size):
            #
            poly = Zx.random_element(euler_phi_pm - 1)  
            
            #
            poly_fast = fast_trace_pm_to_pm_1(p, m, poly)
            poly_calculated = trace_pm_to_pm_1(p, m, poly)   

            #
            assert poly_calculated == poly_fast, "Test {}/{} Failed \n calculated {}  \n fast calculated {}".format(test_number+1, test_size, poly_calculated, poly_fast)
    
        end_time  = perf_counter()
        #print("Test completed, time elapsed: {:2f}".format(end_time - start_time))
    
    def test_multiple(p_range, m_max, tests = 10, test_size = 10):
        #print("\n######\n Testing trace_pm_to_pm_1  \n######\n")
        for i in range(tests):
            p = random_prime(p_range)
            m = ZZ.random_element(1, m_max)
            #print("Teste {}:".format(i+1))
            #print("Prime p = {}, exponent m = {}".format(p, m))
            test(p, m, test_size)
        print("trace pm to pm-1 working as intended")

    test_multiple(p_range = 5, m_max = 4, tests = 5, test_size = 2)

test_pm_to_1()
test_pm_to_pm_1()
