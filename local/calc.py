from ase.calculators.vasp import Vasp2, VaspChargeDensity
from ase.dft.bandgap import bandgap
from ase.io import read
from ase.lattice import BravaisLattice
from matplotlib import pyplot as plt


calc = Vasp2(restart=True, directory="./band/")
bs = calc.band_structure()
# band_structure = BandStructure.read(par_dir + "/band/bs.json")  # 读取上一步得到的能带数据
ref = bs.reference  # 获取费米能级
bs.plot(filename="./bs.jpg", emin=-13, emax=13)  # 画出费米能级
plt.clf()
""":type:BandStructure"""
# band_structure = band_structure.subtract_reference()  # 得到减去费米能级的能带数据
gap = bandgap(calc)
vcd = VaspChargeDensity("./scf/LOCPOT")
local_potential = (vcd.chg[-1] * vcd.atoms[-1].get_volume()).mean(axis=(0, 1))
plt.plot(local_potential)
plt.savefig("./local_potential.jpg")
VBM = local_potential[0] - ref
CBM = VBM + gap[0]
with open("./log.txt", "a+") as tmp:
    print(f"gap:{gap[0]},VBM:{VBM},CBM:{CBM}", file=tmp)
