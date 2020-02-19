# GAMESS RI-MP2 mini-app


General Atomic and Molecular Electronic Structure System (GAMESS) 
is a popular quantum chemistry software package 
which has been around since the 1980s. 
It can calculate a wide variety of molecular properties 
using electronic structure methods. 
One of the methods implemented in GAMESS is 
resolution of identity Moller-Plesset perturbation (RI-MP2) theory. 
RI-MP2 is an electron correlation method, 
which is a class of methods 
that include instantaneous electron-electron interactions, 
and are required to perform accurate energy and property calcula- tions for certain classes of molecular systems. 
Of the electron correlation methods, 
RI-MP2 tends to be one of the more computationally inexpensive methods, 
but the formal computational complexity is still O(N<sup>5</sup>), 
where N is a measure of system size.   


The __GAMESS RI-MP2 mini-app__ computes 
__the correlation energy__
with the Hartree-Fock energy and wave-function given as inputs. 
The inputs were generated from GAMESS. 

## Inputs for GAMESS RI-MP2 mini-app
In this git repository, there is only one input file, benz.kern.
It is the smallest input. You can find bigger inputs at the following link:  
[https://anl.box.com/v/GAMESS-RI-MP2-Inputs](https://anl.box.com/v/GAMESS-RI-MP2-Inputs)    
On NVIDIA V100 GPUs, we recommend to use c60, w30, or w60 inputs to see meaningful speedups. 

## Running GAMESS RI-MP2 mini-app

### SUMMIT at OLCF

#### Build the executables
    rimp2-cublas:   rimp2 with OpenMP offloading + cublas on GPU,
    rimp2-cublasxt: rimp2 with OpenMP offloading + cublasxt on GPU,
    rimp2-nvblas:   rimp2 with OpenMP offloading + nvblas on GPU,
    rimp2-essl:     rimp2 with OpenMP threading  + ESSL on CPU, and
    rimp2-serial:   rimp2 with a single thread   + ESSL on CPU.
    
    $ source source_me_OLCF
    $ make clean
    $ make all
    

#### Run via an interactive job:
    $ bsub -P <your project code> -nnodes 1 -W 120 -Is /bin/bash
    $ source source_me_OLCF
    $ NMPI=x INPUT=xxx EXEC='rimp2-xxx rimp2-yyy' ./run_gpu.sh             
    $ NMPI=x NTHREAD=x INPUT=xxx EXEC='rimp2-zzz' ./run_cpu.sh
        # NMPI is the number of MPIs. If it doesn't exist, NMPI is set to 1.
        # NTHREAD is the number of OpenMP threads per MPI. If it doesn't exist, NTHREAD is set to min(42, 42*NNODES/NMPI).
        # INPUT is the input name. If it doesn't exist, INPUT is set to benz.kern.
        # EXEC is the executable name(s). If it doesn't exist. EXEC is set to 'rimp2-cublas rimp2-cublasxt rimp2-nvblas' for run_gpu.sh, and 'rimp2-essl rimp2-serial' for run_cpu.sh

#### Run via a batch job:
    $ bsub run_batch_OLCF_example.sh            
        # This example runs rimp2_gpu.sh and rimp2_cpu.sh (only with rimp2-essl) with two inputs (cor.kern, and c60.kern)
        #     on 4 SUMMIT nodes with 1, 2, 4, 6, 12, and 24 MPI ranks ( 1 GPU/MPI, 7 CPU threads/MPI).
        # You may modify this example script for your own tests.



### JLSE Skylake nodes at ALCF

#### Build the executables 
    rimp2-mkl:      rimp2 with OpenMP threading + MKL on CPU

    $ source source_me_JLSE_Intel
    $ make clean
    $ make all

#### Run via an interactive job:
    $ qsub -I -n 1 -t 120 -q skylake_8180
    $ source source_me_JLSE_Intel
    $ NMPI=x NTHREAD=x INPUT=xxx EXEC='rimp2-zzz' ./run_cpu.sh
        # NMPI is the number of MPIs. If it doesn't exist, NMPI is set to 1.
        # NTHREAD is the number of OpenMP threads per MPI. If it doesn't exist, NTHREAD is set to min(56, 56*NNODES/NMPI).
        # INPUT is the input name. If it doesn't exist, INPUT is set to benz.kern.
        # EXEC is the executable name(s). If it doesn't exist. EXEC is set to 'rimp2-mkl' for run_cpu.sh

#### Run via a batch job:
    $ qsub ./run_batch_JLSE_example.sh            
        # This example runs rimp2_cpu.sh with two inputs (cor.kern, and c60.kern)
        #     on 1 Skylake 8180 node with 1, 2, and 4 MPI ranks ( 56 threads in total).
        # You may modify this example script for your own tests.


