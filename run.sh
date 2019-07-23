rm -rf Makefile
ln -s Makefile_CPU Makefile
make clean
make
cp rimp2 rimp2-cpu
OMP_PROC_BIND=spread OMP_NUM_THREADS=42 jsrun -n 1 -c 42 -a 1 -b packed:42 -g 0 ./rimp2-cpu $1


rm -rf Makefile
ln -s Makefile_SER Makefile
make clean
make
cp rimp2 rimp2-ser
OMP_PROC_BIND=spread OMP_NUM_THREADS=42 jsrun -n 1 -c 42 -a 1 -b packed:42 -g 0 ./rimp2-ser $1


