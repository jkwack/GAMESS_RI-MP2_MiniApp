//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include "mkl.h"
#include "common.h"
#define QVV(I,J) QVV[I*NVIR+J]

void RIMP2_Energy_Whole_Combined(double *E2){

    double *QVV;
    double E2_local=0.0E0;

    int Nthreads=omp_get_max_threads();
    #pragma omp parallel num_threads(Nthreads) shared(B32,eij,eab,NAUXBASD,NACT,NVIR,E2) private(QVV,E2_local)
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
           dgemm("T", "N",
                &m,     &n,     &k,
                &alpha, &B32[IACT*NAUXBASD*NVIR],  &lda,
                        &B32[JACT*NAUXBASD*NVIR],  &ldb,
                &beta,  QVV,                       &ldc);


           // Accumulate E2
           double E2_t=0.0;
           for(int IB=0; IB<NVIR; IB++){
           for(int IA=0; IA<NVIR; IA++){
              double Tijab = QVV(IB,IA) / ( eij(JACT,IACT) - eab(IB,IA) );
              double Q_t = 2*QVV(IB,IA);
              E2_t += Tijab * (Q_t - QVV(IA,IB) );
           }}   // loop for IA and IB
           double FAC = (IACT==JACT) ? (1.0E0) : (2.0E0);
           E2_local +=  FAC*E2_t;
       }}   // loop for IACT and JACT

       #pragma omp atomic
       *E2 = *E2 + E2_local;
       delete[] QVV;
    } // end of #pragma omp parallel

}

