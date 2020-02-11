#!/bin/bash

tag=$(date | awk '{print $6 $2 $3 $4}')
NNODES=$(cat $COBALT_NODEFILE | wc -l)
NCOREperNODE=56                             # For dueal socket Intel Skylake 8180M processors
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
  EXEC='rimp2-mkl'
  echo "   EXEC is set to $EXEC. For another EXEC, use EXEC='x y' before this job script."
fi

if [ "$NQVV" != "" ]; then
  echo "   NQVV is set to $NQVV."
fi

for EXE in $EXEC; do
  echo -e "\n\n[[[Running $EXE with $NMPI MPI rank(s) and $NTHREAD threads/MPI ...]]]"
  OMP_PROC_BIND=spread OMP_NUM_THREADS=$NTHREAD mpirun -n $NMPI ./$EXE $INPUT $NQVV
  echo -e "\n\n"
done

