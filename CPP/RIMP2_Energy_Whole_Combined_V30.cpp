//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include "common.h"
#define QVV(I,J) QVV[I*NVIR+J]

#include "mkl.h"
#include "mkl_omp_offload.h"

void RIMP2_Energy_Whole_Combined(double *E2){

    double *QVV;
    double E2_local=0.0E0;
    int dnum=0;
    int nQVV=NVIR*NVIR;
    double *A_tmp, *B_tmp;

    #pragma omp target enter data map(alloc:QVV[0:nQVV]) device(dnum)
    #pragma omp target enter data map(to:NAUXBASD,NACT,NVIR,B32size) device(dnum)
    #pragma omp target enter data map(to:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)
    {
       QVV = new double[NVIR*NVIR];
       E2_local = 0.0E0;

       #pragma omp for schedule(dynamic) 
       for(int JACT=0;JACT<NACT;JACT++){
       for(int IACT=0;IACT<JACT;IACT++){

           // Compute QVV
           int m=NVIR;
           int n=NVIR;
           int k=NAUXBASD;
           int lda=NAUXBASD;
           int ldb=NAUXBASD;
           int ldc=NVIR;
           double alpha = 1.0;
           double beta = 0.0;
           A_tmp = &B32[IACT*NAUXBASD*NVIR];
           B_tmp = &B32[JACT*NAUXBASD*NVIR];
           #pragma omp target data map(to:A_tmp[0:m*lda],B_tmp[0:n*ldb])
           #pragma omp target variant dispatch use_device_ptr(A_tmp,B_tmp,B32,QVV) device(dnum)
           dgemm("T", "N",
                &m,     &n,     &k,
                &alpha, A_tmp,  &lda,
                        B_tmp,  &ldb,
                &beta,  QVV,    &ldc);

           // Accumulate E2
           double E2_t=0.0;
           #pragma omp target teams distribute parallel for reduction(+:E2_t) map(tofrom:E2_t) collapse(2) device(dnum)
           {
              for(int IB=0; IB<NVIR; IB++){
              for(int IA=0; IA<NVIR; IA++){
                 double Tijab = QVV(IB,IA) / ( eij(JACT,IACT) - eab(IB,IA) );
                 double Q_t = 2*QVV(IB,IA);
                 E2_t += Tijab * (Q_t - QVV(IA,IB) );
              }}   // loop for IA and IB
           }   // omp target teams distribute parallel for
           double FAC = (IACT==JACT) ? (1.0E0) : (2.0E0);
           E2_local +=  FAC*E2_t;
       }}   // loop for IACT and JACT

       *E2 = *E2 + E2_local;
       delete[] QVV;
    }

    #pragma omp target exit data map(release:QVV[0:nQVV]) device(dnum)
    #pragma omp target exit data map(release:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)

}

