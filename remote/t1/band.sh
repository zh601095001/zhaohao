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
path="struct"
mkdir $path
submit $path
