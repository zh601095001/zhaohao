from ase.build import surface
from ase.io import read

atoms = read("MoS2_mp-1434_primitive.cif")
atoms.center(10,axis=2)

atoms.edit()
