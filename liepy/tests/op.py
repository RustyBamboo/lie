#!python

import liepy as su
import numpy as np

basis = su.gen_gellmann(3)

basis = [-1j * b for b in basis]

order = [2, 1, 4, 0, 7, 6, 3, 5]
order = np.argsort(order)

basis = [basis[i] for i in order]

d_abc = su.get_d_coefficients(basis)
f_ijk = su.get_structure_constants(basis)


a = np.array([[1, 0, 0, 0, 0, 0, 0, 0]], dtype=np.complex128).T
b = np.array([[0, 0, 0, 0, 0, 0, 0, 1]], dtype=np.complex128).T

# 1,5; 1,6; 1,7; 1,8

assert np.allclose(su.dot(a, b, d_abc), su.dot(b, a, d_abc))
assert np.allclose(su.cross(a, b, f_ijk), -su.cross(b, a, f_ijk))

print(su.dot(a,b, d_abc))

print(d_abc)
