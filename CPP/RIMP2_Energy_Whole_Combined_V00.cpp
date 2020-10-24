//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include "common.h"

void RIMP2_Energy_Whole_Combined(double *E2){

    double *QVV;
    double E2_local;
    int nQVV=NVIR*NBLK*NVIR;
    int iQVV;
    double *A_tmp,*B_tmp;
    int m,n,k;
    int lda,ldb,ldc;
    double alpha = 1.0E0;
    double beta = 0.0E0;

    QVV = new double[nQVV];
    E2_local = 0.0E0;

    for(int JACT=0;JACT<NACT;JACT++){
        for(int IACTmod=0;IACTmod<=JACT/NBLK;IACTmod++){

            // Set a length of blocks in a raw according NBLK
            int IACT = IACTmod*NBLK;
            if((IACTmod+1)*NBLK>JACT+1) { iQVV = JACT - (IACTmod)*NBLK + 1; }
            else { iQVV = NBLK; }

            // Compute QVV using dgemm
            m=NVIR*iQVV;
            n=NVIR;
            k=NAUXBASD;
            lda=NAUXBASD;
            ldb=NAUXBASD;
            ldc=NVIR*iQVV;
            A_tmp = &B32[IACT*NAUXBASD*NVIR];
            B_tmp = &B32[JACT*NAUXBASD*NVIR];

            for (int j = 1; j <= n; ++j) {
            for (int i = 1; i <= m; ++i) {
               double temp = 0.;
               for (int l = 1; l <= k; ++l) {
                  temp += B32[IACT*NAUXBASD*NVIR + l-1 + (i-1)*lda]
                         *B32[JACT*NAUXBASD*NVIR + l-1 + (j-1)*ldb];
               }
               QVV[i-1 + (j-1)*ldc] = temp;
            }}


            // Accumulate E2
            for(int IC=0; IC<iQVV; IC++){
                double E2_t=0.0;
                for(int IB=0; IB<NVIR; IB++){
                    for(int IA=0; IA<NVIR; IA++){
                        double Tijab = 
                            QVV[IA+ IC*NVIR+ IB*NVIR*iQVV]
                            / ( eij[IACT+IC + JACT*NACT]
                                - eab[IA + IB*NVIR] );
                        double Q_t = 
                            QVV[IA+ IC*NVIR+ IB*NVIR*iQVV]
                            + QVV[IA+ IC*NVIR+ IB*NVIR*iQVV];
                        E2_t = E2_t + Tijab * (Q_t - QVV[IB+ IC*NVIR+ IA*NVIR*iQVV]);
                    }   // loop for IA
                }   // loop for IB
                double FAC = (IACT==JACT) ? (1.0E0) : (2.0E0);
                E2_local = E2_local + FAC*E2_t;
            }   // loop for IC
        }   // loop for IACTmod
    }   // loop for JACT

    *E2 = *E2 + E2_local;

    delete[] QVV;
}


