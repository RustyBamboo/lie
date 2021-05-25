#!python

import liepy as su

basis = su.gen_gellmann(3)

f_ijk = su.get_structure_constants(basis)
print(f_ijk)

print(f_ijk[(0,1)])
