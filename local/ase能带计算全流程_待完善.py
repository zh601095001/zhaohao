import os
from shutil import rmtree, copytree, copy, move

from ase.calculators.vasp import Vasp2, VaspChargeDensity
from ase.dft.bandgap import bandgap
from ase.io import read
from ase.spectrum.band_structure import BandStructure
from matplotlib import pyplot as plt


class Screen:
    def __init__(self):
        self.file_info = None
        self.command = "mpirun -n 24 /home/share/software/vasp.5.4.1.05Feb16/intel/vtst/vasp_std"

    def build_folders(self, suffix=".cif", folders=("struct", "scf", "band", "calc")):
        path = os.getcwd()
        files = []
        names = []
        file_info = {}
        for file in os.listdir(path):
            if os.path.isfile(file):
                name, _suffix = os.path.splitext(file)
                if _suffix == suffix:  # 判断文件后缀是否符合条件
                    files.append(file)
                    os.mkdir(name)
                    names.append(name)
                    for folder in folders:
                        os.mkdir(f"{name}/{folder}")
        for folder in folders:
            file_info[folder] = {}
        for folder in folders:
            for name in names:
                file_info[folder][name + suffix] = name + "/" + folder
        self.file_info = file_info

    def struct(self):
        struct_file = self.file_info["struct"]
        for filename, path in struct_file.items():
            atoms = read(filename)
            move(filename, os.path.dirname(path))
            calc = Vasp2(command=self.command, directory=path)
            atoms.calc = calc
            calc.set(
                # K-POINTS
                # kpts=(6, 6, 1), gamma=False,
                kspacing=0.5, kgamma=".TRUE.",
                # POTCAR
                xc="hse06",
                # INCAR
                ismear=0, sigma=0.05, prec="Accurate", encut=450.0, lreal="Auto", algo="Fast", ibrion=2,
                nsw=500, ediffg=-0.02, ediff=1E-6,
                lwave=False, lcharg=False, isif=3
            )
            atoms.get_potential_energy()
            atoms.write(path + "/atoms.traj")

    def scf(self):
        scf_files = self.file_info["scf"]
        for filename, path in scf_files.items():
            par_dir = os.path.dirname(path)
            atoms = read(par_dir + "/struct/atoms.traj")
            calc = Vasp2(command=self.command, directory=path)
            atoms.calc = calc
            calc.set(
                # K-POINTS
                # kpts=(6, 6, 1), gamma=False,
                kspacing=0.5, kgamma=".TRUE.",
                # POTCAR
                xc="PBE",
                # INCAR
                ismear=0, sigma=0.05, prec="Accurate", encut=450.0, lreal="Auto", algo="Normal", ibrion=-1, nsw=0,
                ediffg=-0.02, isym=0,
                lwave=False, lcharg=True, laechg=True, lorbit=11, lvhar=True
            )
            atoms.get_potential_energy()
            atoms.write(path + "/atoms2.traj")

    def band(self):
        band_files = self.file_info["band"]
        for filename, path in band_files.items():
            par_dir = os.path.dirname(path)
            self.copy_search_file(par_dir + "/scf", path)
            atoms = read(path + "/atoms2.traj")
            calc2 = Vasp2(command=self.command, directory=path)
            atoms.calc = calc2
            special_path = atoms.cell.get_bravais_lattice().special_path.split(",")[0]
            calc2.set(
                # K-POINTS
                kpts={"path": special_path, "density": 10}, gamma=True,
                # POTCAR
                xc="PBE",
                # INCAR
                ismear=0, sigma=0.05, prec="Accurate", encut=450.0, lreal="Auto", ibrion=-1, nsw=0,
                ediffg=-0.02, isym=0,
                lwave=False, lcharg=True, laechg=True, lorbit=11, lvhar=True, icharg=11,
                # hse06
                lhfcalc=".TRUE.", hfscreen=0.2, lmaxfock=4, algo="Damped", time=0.4, precfock="Normal"
            )
            atoms.get_potential_energy()
            calc2.band_structure().write(path + "/bs.json")

    def calc(self):
        calc_files = self.file_info["calc"]
        for filename, path in calc_files.items():
            par_dir = os.path.dirname(path)
            calc = Vasp2(restart=True, directory=par_dir + "/band/")
            bs = calc.band_structure()
            # band_structure = BandStructure.read(par_dir + "/band/bs.json")  # 读取上一步得到的能带数据
            ref = bs.reference  # 获取费米能级
            bs.plot(filename=path + "/bs.jpg", emin=-13, emax=13, show=False)  # 画出费米能级
            plt.clf()
            """:type:BandStructure"""
            # band_structure = band_structure.subtract_reference()  # 得到减去费米能级的能带数据
            gap = bandgap(calc)
            vcd = VaspChargeDensity(par_dir + "/scf/LOCPOT")
            local_potential = (vcd.chg[-1] * vcd.atoms[-1].get_volume()).mean(axis=(0, 1))
            plt.plot(local_potential)
            plt.savefig(path + "/local_potential.jpg")
            VBM = local_potential[0] - ref
            CBM = VBM + gap[0]
            with open(path + "/log.txt", "a+") as tmp:
                print(f"path:{os.path.dirname(path)},gap:{gap[0]},VBM:{VBM},CBM:{CBM}", file=tmp)

    @staticmethod
    def copy_search_file(src_dir, des_dir):
        ls = os.listdir(src_dir)
        for line in ls:
            file_path = os.path.join(src_dir, line)
            if os.path.isfile(file_path):
                copy(file_path, des_dir)


if __name__ == '__main__':
    test = Screen()
    test.build_folders()
    test.struct()
    test.scf()
    test.band()
    test.calc()
