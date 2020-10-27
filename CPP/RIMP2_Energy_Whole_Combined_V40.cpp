//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include "mkl.h"
#include "mkl_omp_offload.h"
#include "common.h"
#define QVV(I,J) QVV[I*NVIR+J]

void RIMP2_Energy_Whole_Combined(double *E2){

    double *QVV;
    double E2_local=0.0E0;
    int dnum=0;
    QVV = new double[NVIR*NVIR];
    double *B32I, *B32J;

    #pragma omp target enter data map(alloc:QVV[0:NVIR*NVIR]) device(dnum)
    #pragma omp target enter data map(to:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)
    for(int JACT=0;JACT<NACT;JACT++){
    for(int IACT=0;IACT<=JACT;IACT++){

        // Compute QVV
        int n=NVIR;
        int k=NAUXBASD;
        double one = 1.0;
        double zero = 0.0;
        B32I = &B32(IACT,0,0);
        B32J = &B32(JACT,0,0);
//        #pragma omp target data map(to:B32I[0:n*k],B32J[0:n*k])
        #pragma omp target variant dispatch use_device_ptr(B32I,B32J,B32,QVV) device(dnum)
        dgemm("T","N",&n,&n,&k,&one,B32I,&k,B32J,&k,&zero,QVV,&n);

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

    #pragma omp target exit data map(release:QVV[0:NVIR*NVIR]) device(dnum)
    #pragma omp target exit data map(release:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)

}

