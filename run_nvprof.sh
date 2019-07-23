tag=$(date | awk '{print $6 $2 $3 $4}')

for EXE in rimp2-cublasxt rimp2-cublas rimp2-nvblas; do
  echo -e "\n\n[[[Running $EXE with nvprof ...]]]"
  OMP_NUM_THREADS=1 jsrun -n 1 -c 1 -g 1 --smpiargs="none" nvprof -s -f -o nvprof_timeline_${EXE}_${tag}_${1}.nvvp ./${EXE} $1 $2
  echo -e "\n\n"
done



