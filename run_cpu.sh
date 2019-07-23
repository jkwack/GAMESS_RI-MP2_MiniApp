tag=$(date | awk '{print $6 $2 $3 $4}')

for EXE in rimp2-cpu rimp2-ser; do
  echo -e "\n\n[[[Running $EXE ...]]]"
  OMP_PROC_BIND=spread OMP_NUM_THREADS=42 jsrun -n 1 -c 42 -a 1 -b packed:42 -g 0 ./$EXE $1 $2
  echo -e "\n\n"
done



