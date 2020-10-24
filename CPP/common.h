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
extern int NBLK;
extern int B32size;
extern double E2_ref;
#endif
