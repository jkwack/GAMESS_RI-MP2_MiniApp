#!/bin/bash
### Begin BSUB options
#BSUB -P CSC354
#BSUB -J GAMESS-RIMP2-BATCH
#BSUB -W 02:00
#BSUB -nnodes 4
### End BSUB options and begin shell commands

source source_me_OLCF

for iINPUT in cor.kern c60.kern; do
  for iMPI in 24 12 6 4 2 1; do
    NMPI=$iMPI INPUT=$iINPUT ./run_gpu.sh > ${LSB_JOBNAME}.${LSB_JOBID}.GPU.NMPI=${iMPI}.INPUT=${iINPUT}
    NMPI=$iMPI INPUT=$iINPUT NTHREAD=7 EXEC='rimp2-essl' ./run_cpu.sh > ${LSB_JOBNAME}.${LSB_JOBID}.CPU.NMPI=${iMPI}.NTHREAD=7.INPUT=${iINPUT}
  done
done

