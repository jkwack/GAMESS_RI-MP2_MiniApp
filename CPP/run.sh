module load nvhpc

nvc++ -mp=gpu -std=c++11 -DOMP -DOFFLOAD -DCUDA_VERSION -I/soft/compilers/nvhpc/Linux_x86_64/21.7/cuda/11.4/targets/x86_64-linux/include  RIMP2_Energy_Whole_Combined_V45_cublas.cpp GAMESS-RIMP2-CorrEng-Tutorial.c\
pp -lcublas -lcudart -L/soft/compilers/nvhpc/Linux_x86_64/21.7/math_libs/11.4/lib64/ -L/soft/compilers/nvhpc/Linux_x86_64/21.7/cuda/11.4/targets/x86_64-linux/lib

export OMP_TARGET_OFFLOAD=mandatory
export OMP_NUM_THREADS=1

ncu -f --set detailed -k regex:"nvkernel__Z27RIMP2_Energy_Whole_CombinedPd_F1L49_1" -o output_A100_ncu_2 ./a.out w30.rand
