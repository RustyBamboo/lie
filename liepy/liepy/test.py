import algebra

import numpy as np

d = 3

alg = algebra.Algebra_su(d) 

w = np.exp(2j * np.pi/d)
i, j = np.meshgrid(np.arange(d), np.arange(d))
H = np.power(w, i * j ) / np.sqrt(d)

T = np.diag([np.exp(np.pi * 1j * s**3/ d**2) for s in range(d)])

H = algebra.Gate_su(H, alg)
T = algebra.Gate_su(T, alg)

def commu(A, B):
    return A.dot(B) - B.dot(A)

def anti(A,B):
    return A.dot(B) + B.dot(A)

print(np.isclose(H @ T, commu(H.vec(), T.vec())))

print(np.isclose(anti(H.vec(), T.vec()), H * T))

print(H ** T)
