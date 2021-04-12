 ## Run the CPU code

 ```
 $ make my_rimp2_cpu
 $ OMP_PROC_BIND=spread OMP_NUM_THREADS=42 jsrun -n 1 -c 42 -a 1 -b packed:42 -g 0 ./my_rimp2_cpu w30
 ```

 ## Run the CPU serial code with NVBLAS

 ```
 $ make my_rimp2
 $ OMP_NUM_THREADS=1 jsrun -n 1 -c 1 -g 1 ./my_rimp2 w30
 ```

 ## Run the GPU code with NVBLAS

```
$ make my_rimp2_v2
$ OMP_NUM_THREADS=1 jsrun -n 1 -c 1 -g 1 ./my_rimp2_v2 w30
```