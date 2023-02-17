#!/bin/bash

tag=$(date | awk '{print $6 $2 $3 $4}')
NNODES=$SLURM_JOB_NUM_NODES
MAXGPU=$(($NNODES*8))

echo "Running this script with $NNODES node(s) with up to $MAXGPU GPUs:"

if [ "$NMPI" != "" ]; then
  if [ "$NMPI" -gt "$MAXGPU" ]; then
    NMPI=$MAXGPU
    echo "   Too many MPI ranks are requested, so it is adjusted to "$NMPI"."
  fi
  echo "   NMPI is set to $NMPI."
else
  NMPI=1
  echo "   NMPI is set to $NMPI. For another NMPI, use NMPI=x before this job script."
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
  EXEC='rimp2-hip rimp2-serial'
  echo "   EXEC is set to $EXEC. For another EXEC, use EXEC='x y' before this job script."
fi

if [ "$NQVV" != "" ]; then
  echo "   NQVV is set to $NQVV."
fi

for EXE in $EXEC; do
  echo -e "\n\n[[[Running $EXE with $NMPI MPI rank(s)...]]]"
  set -x
  time -p OMP_NUM_THREADS=1 srun -N${NNODES} -n${NMPI} -c1 --gpus-per-node=8 --gpu-bind=closest ./$EXE $INPUT $NQVV
  set +x
  #OMP_NUM_THREADS=1 jsrun -n $NMPI -c 1 -g 1  ./$EXE $INPUT $NQVV
  echo -e "\n\n"
done

