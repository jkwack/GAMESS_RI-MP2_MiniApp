#!/bin/bash

tag=$(date | awk '{print $6 $2 $3 $4}')
NNODES=$((($LSB_DJOB_NUMPROC-1)/42))
MAXGPU=$(($NNODES*6))
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
  EXEC='rimp2-cublasxt rimp2-cublas rimp2-nvblas'
  echo "   EXEC is set to $EXEC. For another EXEC, use EXEC='x y' before this job script."
fi

if [ "$NQVV" != "" ]; then
  echo "   NQVV is set to $NQVV."
fi

OUTDIR=/gpfs/alpine/csc354/proj-shared/jkwack

for EXE in $EXEC; do
  echo -e "\n\n[[[Running $EXE with $NMPI MPI rank(s)...]]]"
  OMP_NUM_THREADS=1 jsrun -n $NMPI -c 1 -g 1 --smpiargs="none" nvprof -s -f -o ${OUTDIR}/nvprof_timeline_${EXE}_${tag}_${INPUT}.nvvp ./$EXE $INPUT $NQVV
  echo -e "\n\n"
done


