//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include <iostream>
#include <string>
#include <cstdlib>
#include <new>
#include <fstream>
#include <ctime>

#if defined(OMP) || defined(OFFLOAD)
#include <omp.h>
inline double timer() {return (omp_get_wtime());}
#else
inline double timer() {return (std::clock()/(double)CLOCKS_PER_SEC);}
#endif

#define TOL 1.0E-6

#ifndef COMMON_H_
#define COMMON_H_
extern double *eij, *eab, *B32;           // They were 2D arrays in Fortran
extern double *EIG;                       // 1D array
extern int NAUXBASD,NCOR,NACT,NVIR,NBF;
extern int B32size;
extern double E2_ref;
#endif

#define eab(I,J) eab[I*NVIR+J]
#define eij(I,J) eij[I*NACT+J]
#define B32(I,J,K) B32[I*NAUXBASD*NVIR + J*NAUXBASD + K]

#ifdef CUDA_VERSION
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
  extern cublasHandle_t handle;
#endif
