#!/bin/bash 
### Begin BSUB Options 
#BSUB -P CHM135
#BSUB -J rimp2-w30_w60-32nodes
#BSUB -W 2:00
#BSUB -nnodes 32
##BSUB -alloc_flags "smt4"
### End BSUB Options and begin shell commands


tag=$(date | awk '{print $6 $2 $3 $4}')
INP_DIR=/gpfs/alpine/chm135/proj-shared/buu/fortranKerns/inputKernels

module load cuda
module load essl/6.2.0-20190419

#EXE=rimp2-cublasxt

for EXE in rimp2-cublasxt rimp2-cublas; do
for INP in w30 w60; do
for nMPI in 6 12 24 48 96 192; do
  echo -e "\n\n[[[Running $EXE $INP w/ nMPI=$nMPI ...]]]"
  OMP_NUM_THREADS=1 jsrun -n $nMPI -a 1 -c 1 -g 1 -r 6  ./$EXE ${INP_DIR}/${INP}.kern
  echo -e "\n\n"
done
done
done



