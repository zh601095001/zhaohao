#!/bin/bash

#PBS -N zh
#PBS -l nodes=1:ppn=28
#PBS -q batch
#PBS -V
# /home/share/software/vasp.5.4.1.05Feb16/intel/clean/vasp_std
submit() {
  local dir
  dir=$1
  if test -z "${dir}"; then
    cd "${PBS_O_WORKDIR}" || exit
  else
    cd "${PBS_O_WORKDIR}"/"${dir}" || exit
  fi
  local NP
  NP=$(cat "$PBS_NODEFILE" | wc -l)
  NN=$(cat "$PBS_NODEFILE" | sort | uniq | tee /tmp/nodes.$$ | wc -l)
  cat "$PBS_NODEFILE" >/tmp/nodefile.$$
  mpirun -genv I_MPI_DEVICE ssm -machinefile /tmp/nodefile.$$ -n "$NP" /home/share/software/vasp.5.4.1.05Feb16/intel/clean/vasp_std
  rm -rf /tmp/nodefile.$$
  rm -rf /tmp/nodes.$$
  cd ..
}

cd "${PBS_O_WORKDIR}" || exit
strpath="struct"
scfpath="scf"
bandpath="band"
mkdir $strpath
mkdir $scfpath
mkdir $bandpath

#
# 1.structure optimization
#

cat > INCAR <<!
  INCAR for structure optimization
  ENCUT = 500
  SIGMA = 0.050000
  ADDGRID = T
  LREAL = A
  NELMIN = 5
  EDIFF = 1.00e-04
  EDIFFG = -2.00e-02
  ALGO = N
  #GGA = PE
  PREC = Accurate
  IBRION = 2
  ISIF = 2
  ISMEAR = 0
  NSW = 500
  LCHARG = .FALSE.
  LWAVE = .FALSE.
  LREAL = Auto
!
cat > POSCAR <<!
 MoS2_mp-2815_primitive
   1.000
    3.1902999878000000    0.0000000000000000    0.0000000000000000
   -1.5951499939000000    2.7628808350999998    0.0000000000000000
    0.0000000000000000    0.0000000000000000   18.1298007964999996
   Mo   S
    1    2
Direct
    0.6666700240000000    0.3333300050000000    0.5000049995000000    Mo1
    0.3333300050000000    0.6666700240000000    0.5863149985000000     S1
    0.3333300050000000    0.6666700240000000    0.4136850015000000     S2

!

#####
cat  > KPOINTS <<!
K-Mesh
0
Gamma
   11   11   1
0.0  0.0  0.0
!
vaspkit -task 103
######

# vaspkit -task 102 -kpr 0.04
mv INCAR KPOINTS POSCAR POTCAR -t $strpath
submit $strpath

#
# 2.scf
#

mv ./$bandpath/CONTCAR ./POSCAR
cat > INCAR <<!
  INCAR for scf
  ENCUT = 500
  SIGMA = 0.050000
  EDIFFG = -2.00e-02
  ALGO = Normal
  GGA = PE
  PREC = Accurate
  IBRION = -1
  ISMEAR = 0
  ISYM = 0
  LORBIT = 11
  NSW = 0
  LAECHG = .TRUE.
  LCHARG = .TRUE.
  LVHAR = .TRUE.
  LWAVE = .FALSE.
  LREAL = Auto
!

#####
cat  > KPOINTS <<!
K-Mesh
0
Gamma
   11   11   1
0.0  0.0  0.0
!
vaspkit -task 103
######

# vaspkit -task 102 -kpr 0.04
mv INCAR KPOINTS POSCAR POTCAR -t $scfpath
submit $scfpath

#
# 3.band
#

cat > INCAR <<!
  INCAR for band structure
  ISTART =  1            (Read existing wavefunction; if there)
  ICHARG =  1            (Non-self-consistent: GGA/LDA band structures)
  LREAL  = F             (Projection operators: automatic)
  ENCUT  =  500          (Cut-off energy for plane wave basis set, in eV)
  PREC   =  Accurate       (Precision level)
  LWAVE  = .TRUE.        (Write WAVECAR or not)
  LCHARG = .TRUE.        (Write CHGCAR or not)
  ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)
  Electronic Relaxation
  ISMEAR =  0            (Gaussian smearing; metals:1)
  SIGMA  =  0.05         (Smearing value in eV; metals:0.2)
 # NELM   =  60           (Max electronic SCF steps)
 # NELMIN =  4            (Min electronic SCF steps)
  EDIFF  =  1E-08       (SCF energy convergence; in eV)
  GGA  =  PE             (PBEsol exchange-correlation)
  LORBIT = 11
  NEDOS = 5000
  Ionic Relaxation
 # NELMIN =  6           (Min electronic SCF steps)
  NSW    =  0            (Max electronic SCF steps)
  IBRION =  -1           (Algorithm: 0-MD; 1-Quasi-New; 2-CG)
  ISIF   =  2            (Stress/relaxation: 2-Ions, 3-Shape/Ions/V, 4-Shape/Ions)
  EDIFFG = -0.001         (Ionic convergence; eV/AA)
  ISYM =  0              (Symmetry: 0=none; 2=GGA; 3=hybrids)

#HSE06 Calculation
  LHFCALC= .TRUE.       (Activate HF)
  AEXX   =  0.25        (25% HF exact exchange, adjusted this value to reproduce experimental band gap)
  HFSCREEN= 0.2         (Switch to screened exchange; e.g. HSE06)
  ALGO   =  Damped         (Electronic Minimisation Algorithm; ALGO=58)
  TIME   =  0.4         (Timestep for IALGO5X)
  PRECFOCK= Normal      (HF FFT grid)
  # NKRED    = 2        (Reduce k-grid-even only, see also NKREDX, NKREDY and NKREDZ)
  # HFLMAX   = 4        (HF cut-off: 4d, 6f)
  # LDIAG    = .TRUE.   (Diagnolise Eigenvalues)
!
mv ./$scfpath/POSCAR ./$scfpath/POTCAR ./$scfpath/WAVECAR -t .
vaspkit -task 302
vaspkit -task 251
mv INCAR KPOINTS POSCAR POTCAR WAVECAR KPATH.in -t $bandpath
submit $bandpath

#
# 4.Post-process Band Structure
#
cd ./$bandpath || exit
echo -e "252\n" | vaspkit
mv ./band.py ./$bandpath
NP=$(cat "$PBS_NODEFILE" | wc -l)
NN=$(cat "$PBS_NODEFILE" | sort | uniq | tee /tmp/nodes.$$ | wc -l)
cat "$PBS_NODEFILE" >/tmp/nodefile.$$
mpirun -genv I_MPI_DEVICE ssm -machinefile /tmp/nodefile.$$ -n "$NP" python band.py
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$
cd ..
#
# 5.calculate cbm and vbm
#
cd ./scf || exit
vacuum=$(echo -e "426\n3\n"|vaspkit|grep "^ Vacuum-Level"|cut -d":" -f2|awk '$1=$1')
cd ../band || exit
gap=$(cat BAND_GAP |grep "Band Gap"|cut -d: -f2|awk '$1=$1')
fermi=$(cat BAND_GAP |grep "Fermi Energy"|cut -d: -f2|awk '$1=$1')
vbm=$(echo $(echo "$vacuum-($fermi)")|bc)
cbm=$(echo $(echo "$vbm-$gap")|bc)
cd .. || exit
echo -e "vbm:-$vbm\ncbm:-$cbm\n"
echo -e "vacuum:$vacuum\ngap:$gap\nfermi:$fermi\nvbm:-$vbm\ncbm:-$cbm\n" > result.txt