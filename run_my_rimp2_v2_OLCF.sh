#!/usr/bin/bash
tag=$(date | awk '{print $6 $2 $3 $4}')
input=$1

for EXE in my_rimp2_v2; do
  echo -e "\n\n[[[Running $EXE ...]]]"
  OMP_NUM_THREADS=1 jsrun -n 1 -c 1 -g 1  ./$EXE $input
  echo -e "\n\n"
done


