load("./GSW/GSW.sage")

n = 2 ** 4 * 3 ** 2
B, q = 2, 2**20
sigma = 2
    
gsw = GSW(n, q, sigma)  

def random_msg_gen(B, n):
    k = ZZ.random_element(0, euler_phi(n)+1)
    k1 = ZZ.random_element(1,3)
    if k == euler_phi(n):
        return 0
    else:
        if k1 == 1:
            return Zx(x^k)
        else:
            return Zx(-x^k)

fn = Zx(cyclotomic_polynomial(n))
Zbx = ZZ.quotient(B)['x'] 
Rb = Zbx.quotient(fn)
#print("initial error ", gsw.get_noise(c, Zx(Rb(m).lift()) ))

d = 30
errors = {}
for i in range(d):
    errors[i] = 0

test = 10
for _ in range(test):
    m = random_msg_gen(B, n)
    c = gsw.enc(m)
    for i in range(d):
        m1 = random_msg_gen(B, n)
        c1 = gsw.enc(m1)
        m = m * m1
        c = gsw.mult(c, c1)
        #print("{}: ".format(i), Zx(Rb(m).lift()))
        #print("{}: ".format(i),gsw.dec(c))
        #errors[i]+=(gsw.get_noise(c, Zx(Rb(m).lift())))
        errors[i]+=(gsw.get_noise(c, m % Zx(cyclotomic_polynomial(n))))
        
for i in range(d):
    print((errors[i]/test).n())