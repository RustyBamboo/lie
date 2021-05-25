#!python

import numpy as np
from .liepy import *
from scipy.linalg import expm, logm

class Algebra_su(object):
    def __init__(self, d):
        self.d = d
        self.su_basis = gen_gellmann(self.d)
        self.su_basis_vec = np.column_stack([np.atleast_2d(x.flatten()).T for x in self.su_basis])

        self.f_ikl = get_structure_constants(self.su_basis)
        self.d_abc = get_d_coefficients(self.su_basis)

    def dim(self):
        return self.d

    def __eq__(self, other):
        return self.dim() == other.dim()

    def __str__(self):
        return f"su({self.d}); basis size: {len(self.su_basis)}"

    def __repr__(self):
        return self.__str__()

class Gate_su(object):
    def __init__(self, U, alg):
        self.alg = alg
        self.U = U / np.linalg.det(U) ** (1/alg.dim())
        self.su = logm(self.U) 

        self.L, err, _, _ = np.linalg.lstsq(self.alg.su_basis_vec, self.su.flatten()[:, np.newaxis], rcond=None)
        assert np.allclose(self.su, sum([a * B for a,B in zip(self.L, self.alg.su_basis)]))

    def mat(self):
        return self.U

    def vec(self):
        return self.su

    def __matmul__(self, other):
        assert self.alg == other.alg
        cross = cross(self.L, other.L, self.alg.f_ikl)
        cross = sum([cross[i] * self.alg.su_basis[i] for i in range(len(cross))])
        return cross

    def __mul__(self, other):
        assert self.alg == other.alg
        dot = dot(self.L, other.L, self.alg.d_abc)
        dot = sum([dot[i] * self.alg.su_basis[i] for i in range(len(dot))])

        dot = -self.L.T.dot(other.L) * 4/self.alg.dim() * np.eye(self.alg.dim()) + 2*dot
        return dot

    def __pow__(self, other):
        dot = dot(self.L, other.L, self.alg.d_abc)
        cross = cross(self.L, other.L, self.alg.f_ikl)

        prod = sum([(dot[i] + cross[i]/2) * self.alg.su_basis[i] for i in range(len(dot))])
        prod = -self.L.T.dot(other.L) * 2/self.alg.dim() * np.eye(self.alg.dim())
        return prod

    def __str__(self):
        return f"{self.mat().__str__()}\nIn basis:\n{self.L}"

    def __repr__(self):
        return self.__str__()
        



