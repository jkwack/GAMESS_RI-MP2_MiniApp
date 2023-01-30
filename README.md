
 ## Put the correct nvblas in the path

 ```
 $ source source_me_ALCF_POLARIS
 $ make -f Makefile_polaris clean
 ```

 ## Run the CPU code

 ```
 $ module switch PrgEnv-nvhpc PrgEnv-cray
 $ module unload craype-accel-nvidia80
 $ make -f Makefile_polaris my_rimp2_cpu
 $ OMP_PROC_BIND=spread OMP_NUM_THREADS=64 ./my_rimp2_cpu w30
 ```

 ## Run the CPU serial code with NVBLAS

 ```
 $ module restore
 $ make -f Makefile_polaris my_rimp2
 $ CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=1 ./my_rimp2 w30
 ```

 ## Run the GPU code with NVBLAS

```
 $ module restore
 $ make -f Makefile_polaris my_rimp2_v2
 $ CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=1 ./my_rimp2_v2 w30
```