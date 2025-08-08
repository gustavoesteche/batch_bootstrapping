R = PolynomialRing(ZZ, 'x')
R1 = PolynomialRing(ZZ, 'x1')
R2 = PolynomialRing(ZZ, 'x2')

# experimentando a ida

x = R.gen()
x1 = R1.gen()
x2 = R2.gen()

m1 = 3 ** 1 
m2 = 5 ** 1
m = m1 * m2

f = R(cyclotomic_polynomial(m))

f1 = R1(cyclotomic_polynomial(m1))
poly1 = R1(x1 + 1) % f1

f2 = R2(cyclotomic_polynomial(m2))
poly2 = R2(x2 ** 3 + x2 ** 2 + x2 + 1) % f2
print(poly1)
print(poly2)

i_poly1 = R(poly1(x1^(m/m1)))
i_poly2 = R(poly2(x2^(m/m2))) 
print(i_poly1)
print(i_poly2)

print((R(i_poly1) * R(i_poly2)))


# experimentando a volta

Zx = PolynomialRing(ZZ, 'x')

# definindo o polinomio maior 
m = 15
base = cyclotomic_polynomial(m)
poly = (Zx.random_element(degree = base.degree()) % base)

# Encontrando todos os aneis em que ele vai ser repartido
divisores = list(factor(m))
aneis_polinomiais = []

for prime in divisores:
    aneis_polinomiais.append(cyclotomic_polynomial(prime[0] ** prime[1]))