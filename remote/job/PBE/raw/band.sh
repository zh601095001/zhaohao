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
}

cd "${PBS_O_WORKDIR}" || exit
# structure optimization
strpath="struct"
scfpath="scf"
bandpath="band"
mkdir $strpath
cat > INCAR <<!
  INCAR created by Atomic Simulation Environment
  ENCUT = 450.000000
  SIGMA = 0.050000
  EDIFF = 1.00e-06
  EDIFFG = -2.00e-02
  ALGO = Fast
  GGA = PE
  PREC = Accurate
  IBRION = 2
  ISIF = 3
  ISMEAR = 0
  NSW = 500
  LCHARG = .FALSE.
  LWAVE = .FALSE.
  LREAL = Auto
!
cat > POSCAR <<!
   S Mo
   1.0000000000000000
       6.3659999999999997    0.0000000000000000    0.0000000000000000
      -3.1829999999999985    5.5131177204917368    0.0000000000000000
       0.0000000000000000    0.0000000000000000   15.0000000000000000
    S Mo
     8   4
  Cartesian
    1.5914681700000002  0.9188713304743579  5.9733000000000001
    1.5914681700000002  0.9188713304743579  9.1006499999999999
    4.7744681700000005  0.9188713304743579  5.9733000000000001
    4.7744681700000005  0.9188713304743579  9.1006499999999999
   -0.0000318299999988  3.6754301907202263  5.9733000000000001
   -0.0000318299999988  3.6754301907202263  9.1006499999999999
    3.1829681700000010  3.6754301907202263  5.9733000000000001
    3.1829681700000010  3.6754301907202263  9.1006499999999999
    0.0000318300000004  1.8376875297715107  7.5370499999999998
    3.1830318300000000  1.8376875297715107  7.5370499999999998
   -1.5914681699999986  4.5942463900173793  7.5370499999999998
    1.5915318300000010  4.5942463900173793  7.5370499999999998
!
vaspkit -task 102 -kpr 0.04
mv INCAR KPOINTS POSCAR POTCAR -t $strpath
submit $strpath
# scf
mv CONTCAR ../POSCAR
cd ..
mkdir $scfpath
cat > INCAR <<!
  INCAR created by Atomic Simulation Environment
  ENCUT = 450.000000
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
vaspkit -task 102 -kpr 0.04
mv INCAR KPOINTS POSCAR POTCAR -t $scfpath
submit $scfpath

# band
mv CHGCAR POSCAR POTCAR -t ..
cd ..
mkdir $bandpath
cat > INCAR <<!
  INCAR created by Atomic Simulation Environment
  ENCUT = 450.000000
  SIGMA = 0.050000
  EDIFFG = -2.00e-02
  ALGO = Normal
  GGA = PE
  PREC = Accurate
  IBRION = -1
  ICHARG = 11
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
vaspkit -task 302
cp KPATH.in KPOINTS
mv INCAR KPOINTS POSCAR POTCAR CHGCAR -t $bandpath
submit $bandpath

vaspkit -task 211