#!/bin/bash 
### Begin BSUB Options 
#BSUB -P CHM135
#BSUB -J rimp2-w30_w60-1node-cpu
#BSUB -W 2:00
#BSUB -nnodes 1
##BSUB -alloc_flags "smt4"
### End BSUB Options and begin shell commands


tag=$(date | awk '{print $6 $2 $3 $4}')
INP_DIR=/gpfs/alpine/chm135/proj-shared/buu/fortranKerns/inputKernels

module load cuda
module load essl/6.2.0-20190419


for EXE in rimp2-cpu; do
for INP in w60; do
for nMPI in 4; do
  echo -e "\n\n[[[Running $EXE $INP w/ nMPI=$nMPI ...]]]"
  OMP_PROC_BIND=spread OMP_NUM_THREADS=7 jsrun -n $nMPI -a 1 -c 7 -b packed:7 -g 0 -r $nMPI  ./$EXE ${INP_DIR}/${INP}.kern
  echo -e "\n\n"
done
done
done





