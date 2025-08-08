load("./dual/dbasis.sage")


tests = 10

for _ in range(tests):
    p, m = random_prime(2, 20), ZZ.random_element(1,5)
    #print(p, m)
    #print(dual_factor(p,m))
    basis = [Zx(x^i) for i in range(euler_phi(p**m))]
    dbasis = canon_dbasis(p, m)
    #print(dbasis, basis)
    # deb_canon_dbasis(p, m, dbasis, basis)
    check_dual_basis(basis, dbasis, p, m)

print("The dual basis of the cannonical basis is working as intended")