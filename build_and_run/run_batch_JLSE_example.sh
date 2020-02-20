#!/bin/bash
#COBALT --jobname GAMESS-RIMP2-BATCH
#COBALT -t 02:00
#COBALT -n 1
#COBALT -q skylake_8180
### End COBALT options and begin shell commands

source source_me_JLSE_Intel

NNODES=$(cat $COBALT_NODEFILE | wc -l)
for iINPUT in cor.kern c60.kern; do
  for iMPI in 4 2 1; do
    iTHREAD=$(($NNODES*56/$iMPI))
    NMPI=$iMPI NTHREAD=$iTHREAD INPUT=$iINPUT ./run_cpu.sh > GAMESS-RIMP2-BATCH.${COBALT_JOBID}.CPU.NMPI=${iMPI}.NTHREAD=${iTHREAD}.INPUT=${iINPUT}
  done
done

