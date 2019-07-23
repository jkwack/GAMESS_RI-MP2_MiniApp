for EXE in rimp2-cublasxt rimp2-cublas rimp2-nvblas rimp2-cpu rimp2-ser rimp2-mkl-cpu; do
  echo -e "\n\n\n[[[Building $EXE ...]]]"
  rm -rf $EXE
  rm -rf Makefile
  ln -s Makefile_$EXE Makefile
  make clean
  make
  mv rimp2 $EXE
done




