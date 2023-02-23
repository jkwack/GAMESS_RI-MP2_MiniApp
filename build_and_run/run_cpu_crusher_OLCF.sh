#!/bin/bash

tag=$(date | awk '{print $6 $2 $3 $4}')
NCOREperNODE=128   # For dual-socket AMD EPYC 7A53 64-Core Processor on crusher 
NNODES=$SLURM_JOB_NUM_NODES
MAXTHREAD=$(($NNODES*$NCOREperNODE))
echo "Running this script with $NNODES node(s) with up to $MAXTHREAD CPU threads in total:"

if [ "$NMPI" != "" ]; then
  echo "   NMPI is set to $NMPI."
else
  NMPI=1
  echo "   NMPI is set to $NMPI. For another NMPI, use NMPI=x before this job script."
fi

if [ "$NTHREAD" != "" ]; then
  ALLTHREAD=$(($NMPI*$NTHREAD))
  if [ "$ALLTHREAD" -gt "$MAXTHREAD" ] || [ "$NTHREAD" -gt "$NCOREperNODE" ]; then
    NTHREAD=$(($MAXTHREAD/$NMPI))
    NTHREAD=$(($NTHREAD<=$NCOREperNODE ? $NTHREAD : $NCOREperNODE))
    echo "   Too many CPU threads per rank are requested, so it is adjusted to "$NTHREAD"."
  fi
  echo "   NTHREAD is set to $NTHREAD."
else
  NTHREAD=$(($MAXTHREAD/$NMPI))
  NTHREAD=$(($NTHREAD<=$NCOREperNODE ? $NTHREAD : $NCOREperNODE))
  echo "   NTHREAD is set to $NTHREAD. For another NTHREAD, use NTHREAD=x before this job script."
fi

if [ "$INPUT" == "" ]; then
  INPUT="benz.kern"
  echo "   INPUT is set to $INPUT. For another INPUT, use INPUT=xxxx before this job script."
else
  echo "   INPUT is set to $INPUT."
fi

if [ "$EXEC" != "" ]; then
  echo "   EXEC is set to $EXEC."
else
  EXEC='rimp2-cpu'
  echo "   EXEC is set to $EXEC. For another EXEC, use EXEC='x y' before this job script."
fi

if [ "$NQVV" != "" ]; then
  echo "   NQVV is set to $NQVV."
fi

for EXE in $EXEC; do
  echo -e "\n\n[[[Running $EXE with $NMPI MPI rank(s) and $NTHREAD threads/MPI ...]]]"
  set -x
  OMP_NUM_THREADS=8 srun --distribution=*:block -N${NNODES} -n${NMPI} -c8 ./hello_jobstep
  time -p OMP_NUM_THREADS=8 srun --distribution=*:block -N${NNODES} -n${NMPI} -c8 ./$EXE $INPUT $NQVV
  set +x
  echo -e "\n\n"
done

