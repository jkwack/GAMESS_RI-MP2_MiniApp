//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include "mkl.h"
#include "mkl_omp_offload.h"
#include "common.h"
#define QVV(I,J,K) QVV[I*NVIR*(JACT+1)+J*NVIR+K]

void RIMP2_Energy_Whole_Combined(double *E2){

    double *QVV;
    double E2_local=0.0E0;
    int dnum=0;
    QVV = new double[NVIR*NACT*NVIR];
    double *B32J;
    double opm_time_final = 0;
    #pragma omp target enter data map(alloc:QVV[0:NVIR*NACT*NVIR]) device(dnum)
    #pragma omp target enter data map(to:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)
    for(int JACT=0;JACT<NACT;JACT++){

        // Compute QVV
        const long long  m=NVIR*(JACT+1);
        const long long  n=NVIR;
        const long long  k=NAUXBASD;
        double one = 1.0;
        double zero = 0.0;
        B32J = &B32(JACT,0,0);
        #pragma omp dispatch device(dnum)
        //cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m,n,k,one,B32,k,B32J,k,zero,QVV,m);
	dgemm("T","N",&m,&n,&k,&one,B32,&k,B32J,&k,&zero,QVV,&m);
	double opm_time_init=omp_get_wtime();
        // Accumulate E2
        #pragma omp target teams distribute reduction(+:E2_local) device(dnum)
        for(int IACT=0; IACT<=JACT; IACT++){
           double E2_t=0.0;
           #pragma omp parallel for reduction(+:E2_t)
           for(int IB=0; IB<NVIR; IB++){
           for(int IA=0; IA<NVIR; IA++){
              double Tijab = QVV(IB,IACT,IA) / ( eij(JACT,IACT) - eab(IB,IA) );
              double Q_t = 2*QVV(IB,IACT,IA);
              E2_t += Tijab * (Q_t - QVV(IA,IACT,IB) );
           }}   // loop for IA and IB
           double FAC = (IACT==JACT) ? (1.0E0) : (2.0E0);
           E2_local +=  FAC*E2_t;
        }   // loop for IACT
		opm_time_final=opm_time_final+(omp_get_wtime() - opm_time_init);
    }   // loop for JACT
    printf( "time: %lf\n", opm_time_final);
    *E2 = *E2 + E2_local;
    delete[] QVV;

    #pragma omp target exit data map(release:QVV[0:NVIR*NACT*NVIR]) device(dnum)
    #pragma omp target exit data map(release:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)

}

