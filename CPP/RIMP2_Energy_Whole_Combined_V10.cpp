//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

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
           std::fill_n(QVV,NVIR*NVIR,0.0);
           for (int j = 0; j < NVIR; ++j) {
           for (int i = 0; i < NVIR; ++i) {
           for (int l = 0; l < NAUXBASD; ++l) {
              QVV(j,i) += B32(IACT,i,l)*B32(JACT,i,l);
           }}}


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

