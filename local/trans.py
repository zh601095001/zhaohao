from ase.io import read
from ase.build import bulk
atoms = read("MoS2_mp-1434_primitive.cif")
atoms.center(15,2)
atoms.write("MoS2.cif")