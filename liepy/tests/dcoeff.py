#!python
import numpy as np
import liepy as su
from itertools import permutations

basis = su.gen_gellmann(3)

# Make the matrices match up Wikipedia Gellmann matrices, where T = \lambda / 2
basis = [-1j * b for b in basis]

order = [2, 1, 4, 0, 7, 6, 3, 5]
order = np.argsort(order)

basis = [basis[i] for i in order]

# for b in basis:
#     print(b)

# What the d coeff should be for su(3)
# https://en.wikipedia.org/wiki/Special_unitary_group#Lie_algebra_3
d_abc_exact = {(i, i, 7): 1/np.sqrt(3) for i in range(3)} | {(7,7,7): -1/np.sqrt(3)}
d_abc_exact |= {(i,i,7): -1/(2 * np.sqrt(3)) for i in range(3,8)}
d_abc_exact |= {(2, i, i): 0.5 for i in range(3,5)} | {(2, i, i): -0.5 for i in range(5,7)}
d_abc_exact |= {(1,3,6): -0.5, (0,3,5): 0.5, (0,4,6): 0.5, (1,4,5): 0.5}

d_abc = su.get_d_coefficients(basis)

for (k, v) in sorted(d_abc.items()):
    perm = list(permutations(k))
    found = [p in d_abc_exact for p in perm]
    if any(found):
        i = found.index(True)
        if not np.isclose(v, d_abc_exact[perm[i]]):
            print(f"For {k} expected {d_abc_exact[perm[i]]}, got {v}")
        continue
    else:
        assert False, "Expect d coefficients don't match"

def anti(A,B):
    return A.dot(B) + B.dot(A)

def anti_su(a,b):
    def delta(i,j):
        return 1 if i == j else 0
    I = np.eye(len(basis[a]))

    D = np.zeros_like(I, dtype=np.complex128)
    
    for c in range(8):
        perm = list(permutations((a,b,c)))
        found = [p in d_abc for p in perm]
        if any(found):
            i = found.index(True)
            
            D += 2*d_abc[perm[i]] * basis[c]

    return delta(a,b) * I*4/3 + D


for a in range(8):
    for b in range(8):
        expected = anti(basis[a], basis[b])
        got = anti_su(a,b)
        if not np.allclose(expected, got):
            print(a, b)
            print(expected)
            print(got)
            print("====")
