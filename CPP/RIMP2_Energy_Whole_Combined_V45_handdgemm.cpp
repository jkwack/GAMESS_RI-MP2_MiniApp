//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#ifdef MKL
#include "mkl.h"
#include "mkl_omp_offload.h"
#endif
#include "common.h"
#define QVV(I,J,K) QVV[I*NVIR*(JACT+1)+J*NVIR+K]

void dgemm_hand( int m, int n, int k, double *A, double *B, double *C);

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
        int m=NVIR*(JACT+1);
        int n=NVIR;
        int k=NAUXBASD;
        double one = 1.0;
        double zero = 0.0;

	B32J = &B32(JACT,0,0);

	printf( "%d out of %d\n", JACT, NACT);
	dgemm_hand( m,n,k, B32, B32J, QVV);
	
	// int cublas_error = cublasDgemm(handle,CUBLAS_OP_T, CUBLAS_OP_N, m,n,k, &one, B32, k, B32J, k, &zero, QVV, m);
	// if( cublas_error != CUBLAS_STATUS_SUCCESS )
	//   {
	//     printf( "failed %d.\n", cublas_error );
	//     exit(1);
	//   }
	// cudaDeviceSynchronize();

	//	double opm_time_init=omp_get_wtime();
#pragma omp target update to(QVV[0:NVIR*NACT*NVIR])
        // Accumulate E2
	double opm_time_init=omp_get_wtime();
#pragma omp target teams distribute reduction(+:E2_local)
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

void dgemm_hand( int m, int n, int k, double * __restrict__ A, double * __restrict__ B , double * __restrict__ C)
{

  #pragma omp parallel for
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      //      C[i][j] = 0;
      //      C[i*n+j] = 0;
      C[j*m+i] = 0;
      for (int l = 0; l < k; ++l) {
	//	C[i*n +j] += A[i*k+l]*B[j*k+l];
	C[j*m +i] += A[i*k+l]*B[j*k+l];
	//	C[i][j] += A[i][l]*B[l][j];
	//	C[i][j] += A[i][l]*B[j][l];

      }}}
  
}
