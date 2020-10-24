//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include "common.h"

#if defined(MKL)
#include "mkl.h"
   #if defined(OFFLOAD)
   #include "mkl_omp_offload.h"
   #endif
#endif



#if defined(OFFLOAD)
int dnum=0;
#endif

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

#if defined(OMP)
    int Nthreads=omp_get_max_threads();
    #pragma omp parallel num_threads(Nthreads) shared(NAUXBASD,NCOR,NACT,NVIR,NBF,NBLK,B32size,E2,nQVV,alpha,beta) private(QVV,E2_local,iQVV,A_tmp,B_tmp,m,n,k,lda,ldb,ldc)
    {
#endif
    QVV = new double[nQVV];
    E2_local = 0.0E0;

#if defined(OFFLOAD)
    #pragma omp target enter data map(alloc:QVV[0:nQVV]) device(dnum)
    #pragma omp target enter data map(to:NAUXBASD,NCOR,NACT,NVIR,NBF,NBLK,B32size) device(dnum)
    #pragma omp target enter data map(to:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)
#endif

#if defined(OMP)
    #pragma omp for schedule(dynamic)
#endif
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

#if defined(OFFLOAD)
            #pragma omp target data map(to:A_tmp[0:m*lda],B_tmp[0:n*ldb])
#endif
#if defined(MKL)
   #if defined(OFFLOAD)
            #pragma omp target variant dispatch use_device_ptr(A_tmp,B_tmp,B32,QVV) device(dnum)
   #endif  // for #if defined(OFFLOAD)
            dgemm("T", "N",
                &m,     &n,     &k,
                &alpha, A_tmp,  &lda,
                        B_tmp,  &ldb,
                &beta,  QVV,    &ldc);

#else   // for a non-MKL version
#if defined(V1)
    #if defined(OFFLOAD)
            #pragma omp target teams distribute parallel for collapse(2) device(dnum)
    #endif
            for (int j = 1; j <= n; ++j) {
            for (int i = 1; i <= m; ++i) {
               double temp = 0.;
               for (int l = 1; l <= k; ++l) {
                  temp += B32[IACT*NAUXBASD*NVIR + l-1 + (i-1)*lda]
                         *B32[JACT*NAUXBASD*NVIR + l-1 + (j-1)*ldb];
               }
               QVV[i-1 + (j-1)*ldc] = temp;
            }}
#elif defined(V2_TILE) // blocking
#define TILE 1024
    #if defined(OFFLOAD)
            #pragma omp target teams distribute parallel for collapse(2) device(dnum)
    #endif
            for (int j = 1; j <= n; ++j) {
            for (int i = 1; i <= m; ++i) {
               QVV[i-1 + (j-1)*ldc] = 0.0;
            }}

            for (int jb=1; jb <= n; jb+=TILE) { 
            for (int ib=1; ib <= m; ib+=TILE) {
            for (int lb=1; lb <= k; lb+=TILE) {
               int je=std::min(jb+TILE-1,n);
               int ie=std::min(ib+TILE-1,m);
               int le=std::min(lb+TILE-1,k);
    #if defined(OFFLOAD)
               #pragma omp target teams distribute parallel for simd collapse(2) device(dnum)
    #endif
               for (int j = jb; j <= je; ++j) {
               for (int i = ib; i <= ie; ++i) {
                  double temp = 0.;
                  for (int l = lb; l <= le; ++l) {
                     temp += B32[IACT*NAUXBASD*NVIR + l-1 + (i-1)*lda]
                            *B32[JACT*NAUXBASD*NVIR + l-1 + (j-1)*ldb];
                  }
                  QVV[i-1 + (j-1)*ldc] += temp;
               }}
            }}}
#else
            std::cout<<"Error! Please build this code again with either -DV1 or -DV2_TILE\n";
#endif  // for #if defined(V1)

#endif  // for #if defined(MKL)

//std::cout<<"E2 Accumulation \n";

#if defined(OFFLOAD)
            #pragma omp target map(tofrom:E2_local) device(dnum)
            {
            #pragma omp teams distribute reduction(+:E2_local)
            {
#endif
            // Accumulate E2
            for(int IC=0; IC<iQVV; IC++){
                double E2_t=0.0;
#if defined(OFFLOAD)
                #pragma omp parallel for reduction(+:E2_t) collapse(2)
                {
#endif
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
#if defined(OFFLOAD)
                }   // omp parallel for
#endif    
                double FAC=2.0E0;
                if(IACT+IC == JACT) FAC=1.0E0;
                E2_local = E2_local + FAC*E2_t;
            }   // loop for IC
#if defined(OFFLOAD)
            }   // omp teams distribute
            }   // omp target
#endif
        }   // loop for IACTmod
    }   // loop for JACT

#if defined(OFFLOAD)
    #pragma omp target exit data map(release:QVV[0:nQVV]) device(dnum)
    #pragma omp target exit data map(release:eij[0:NACT*NACT],eab[0:NVIR*NVIR],B32[0:B32size]) device(dnum)
#endif


#if defined(OMP)
    #pragma omp atomic
#endif
    *E2 = *E2 + E2_local;

    delete[] QVV;
#if defined(OMP)
    }    // end of #pragma omp parallel
#endif
}


