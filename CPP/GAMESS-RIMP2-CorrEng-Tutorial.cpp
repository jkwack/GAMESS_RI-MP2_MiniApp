//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include <iostream>
#include <string>
#include <cstdlib>
#include <new>
#include <fstream>
#include <ctime>
#include "mkl.h"
#if defined(OMP) || defined(OFFLOAD)
#include <omp.h>
#endif
#if defined(OFFLOAD)
#include "mkl_omp_offload.h"
int dnum=0;
#endif

// struct rimp2_input {
//     double *eij, *eab;//, *B32;           // They were 2D arrays in Fortran
//     double *EIG;                       // 1D array
//     int NAUXBASD,NCOR,NACT,NVIR,NBF;
//     int NBLK=0;
// //    int B32size;
//     double E2_ref;
// } my;
double *eij, *eab, *B32;           // They were 2D arrays in Fortran
double *EIG;                       // 1D array
int NAUXBASD,NCOR,NACT,NVIR,NBF;
int NBLK=0;
int B32size;
double E2_ref;

void RIMP2_Energy_Whole_Combined(double *E2);
void Initialization(int argc, char *argv[]);
void Finalization();
void Read_Input_File(std::string fname);

#if defined(OMP) || defined(OFFLOAD)
inline double timer() {return (omp_get_wtime());}
#else
inline double timer() {return (std::clock()/(double)CLOCKS_PER_SEC);}
#endif

int main(int argc, char *argv[]){

    // Energy
    double E2, E2_diff, Rel_E2_error;

    // Wall time
    double dt;
    double tic,toc;

#if defined(OMP)
    std::cout<<"You are running the code with OpenMP threading.\n";
#elif defined(OFFLOAD)
    std::cout<<"You are running the code with OpenMP offloading on GPU.\n";
#else
    std::cout<<"You are running the code serially.\n";
#endif

    // Read or generate iput data
    Initialization(argc, argv);

    // Warming up
    E2 = 0.0;
    // Correlation Energy Accumulation
//    RIMP2_Energy_Whole_Combined(E2);

    // Measuing the performance of Correlation Energy Accumulation
    E2 = 0.0;
    tic = timer();
    RIMP2_Energy_Whole_Combined(&E2);
    toc = timer();
    dt = toc-tic;

    // Report the performance data and pass/fail status
    E2_diff = E2 - E2_ref;
    Rel_E2_error = abs(E2_diff/E2_ref);
    std::cout<<"\tResults:\n";
#if defined(OMP)
    std::cout<<"\t\tNumber of OMP threads                   = "<<omp_get_max_threads()<<"\n";
#endif
//    std::cout<<"\t\tReference MP2 corr. energy              = "<<E2_ref<<"\n";
    std::cout<<"\t\tRel. error of computed MP2 corr. energy = "<<Rel_E2_error<<"\n";
    std::cout<<"\t\tWall time                               = "<<dt<<" sec\n";
    if (Rel_E2_error <= 1.0E-6) {
        std::cout<<"\t\tPassed :-) \n";
    }
    else {
        std::cout<<"\t\tFailed :-) \n";
    }

    // Finalize
    Finalization();
    return(0);

}



void RIMP2_Energy_Whole_Combined(double *E2){

    double *QVV;
    double E2_local;
    int nQVV=NVIR*NBLK*NVIR;
    int iQVV;

#if defined(OMP)
    int Nthreads=omp_get_max_threads();
    #pragma omp parallel num_threads(Nthreads) shared(NAUXBASD,NCOR,NACT,NVIR,NBF,NBLK,B32size,E2,nQVV) private(QVV,E2_local,iQVV)
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
            int m=NVIR*iQVV;
            int n=NVIR;
            int k=NAUXBASD;
            double alpha = 1.0E0;
            double beta = 0.0E0;
            int lda=NAUXBASD;
            int ldb=NAUXBASD;
            int ldc=NVIR*iQVV;
#if defined(OFFLOAD)
            #pragma omp target variant dispatch use_device_ptr(B32,QVV) device(dnum)
            {
#endif
            dgemm("T", "N",
                &m,     &n,     &k,
                &alpha, &B32[IACT*NAUXBASD*NVIR], &lda,
                        &B32[JACT*NAUXBASD*NVIR], &ldb,
                &beta,  QVV,    &ldc);
#if defined(OFFLOAD)
            }
#endif

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
std::cout<<E2_local<<"\n";
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


void Initialization(int argc, char *argv[]){

    std::string fname;

    // Read input filename
    if (argc > 1) fname=argv[1];
    if (fname=="benz.kern" || fname=="cor.kern" || fname=="c60.kern" || fname=="w30.kern" || fname=="w60.kern") {
        std::cout<<"\tReading data from "<<fname<<std::endl;
        Read_Input_File(fname);
    }
    else{
        if(fname=="benz.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of benz.kern\n";
            NAUXBASD=420;
            NCOR=6;
            NACT=15;
            NVIR=93;
            NBF=120;
        }
        else if(fname=="cor.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of cor.kern\n";
            NAUXBASD=1512;
            NCOR=24;
            NACT=54;
            NVIR=282;
            NBF=384;
        }
        else if (fname=="c60.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of c60.kern\n";
            NAUXBASD=3960;
            NCOR=60;
            NACT=120;
            NVIR=360;
            NBF=540;
        }
        else if (fname=="w30.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of w30.kern\n";
            NAUXBASD=2520;
            NCOR=30;
            NACT=120;
            NVIR=570;
            NBF=750;
        }
        else if (fname=="w60.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of w60.kern\n";
            NAUXBASD=5040;
            NCOR=60;
            NACT=240;
            NVIR=1140;
            NBF=1500;
        }
        else{
            std::cout<<"\tError!\n\tOne of the followings should be used as an input:\n";
            std::cout<<"\t\tbenz.kern, cor.kern, c60.kern, w30.kern, or w60.kern for actual data sets, or\n";
            std::cout<<"\t\tbenz.rand, cor.rand, c60.rand, w30.rand, or w60.rand for arbitrary data sets.\n";
        }

        // Generate MO energy
        EIG= new double[NBF];
        std::fill_n(EIG,NBF,1);
        for (int ii=0;ii<NACT;ii++) {
            EIG[ii+NCOR] = 2.0;
        }

        // Generate B32
        B32size=NAUXBASD*NVIR*NACT;
        B32= new double[B32size];
        std::fill_n(B32,B32size,1.0/NAUXBASD);

        // Compute the corresponding mp2 corr energy
        E2_ref=0.5*(NVIR*NVIR)*(NACT*NACT)/(NAUXBASD*NAUXBASD);
    }

    // Some parameters
    int NOCC = NCOR + NACT;

    // Virt-Virt MO energy pairs
    eab = new double[NVIR*NVIR];
    for (int IB=0;IB<NVIR;IB++){
        for (int IA=0;IA<=IB;IA++){
            eab[IA+IB*NVIR] = EIG[IA+NOCC] + EIG[IB+NOCC];
            eab[IB+IA*NVIR] = eab[IA+IB*NVIR];
        }
    }

    // occ-occ MO energy pairs
    eij = new double[NACT*NACT];
    for (int JJ=0;JJ<NACT;JJ++){
        for (int II=0;II<=JJ;II++){
            eij[II+JJ*NACT] = EIG[II+NCOR] + EIG[JJ+NCOR];
        }
    }

    // Read the second command line argument for NBLK
    if (argc > 2) {NBLK=std::atoi(argv[2]);}
    else {NBLK=NACT;}

    // Print out the summary of the input
    std::cout<<"\tNAUXBASD NCOR NACT NVIR NBF = "<<NAUXBASD<<" "<<NCOR<<" "<<NACT<<" "<<NVIR<<" "<<NBF<<"\n";
    std::cout<<"\tNBLK = "<<NBLK<<"\n";
    std::cout<<"\tMemory Footprint:\n";
    std::cout<<"\t\tB32[ "<<NAUXBASD*NVIR<<" , "<<NACT<<" ] = "<<NAUXBASD*NVIR*NACT*8.E-6<<" MB\n";
    std::cout<<"\t\teij[ "<<NACT<<" , "<<NACT<<" ] = "<<NACT*NACT*8.E-6<<" MB\n";
    std::cout<<"\t\teab[ "<<NVIR<<" , "<<NVIR<<" ] = "<<NVIR*NVIR*8.E-6<<" MB\n";
    std::cout<<"\t\tQVV[ "<<NVIR<<" , "<<NACT<<" , "<<NVIR<<" ] = "<<NVIR*NACT*NVIR*8.E-6<<" MB\n";

}



void Finalization(){
    delete[] B32;
    delete[] EIG;
    delete[] eij;
    delete[] eab;
}



void Read_Input_File(std::string fname){
    std::ifstream myfile;
    int a;
    double d;

    myfile.open(fname,std::ios::in | std::ios::binary);

    if (myfile.is_open()){

        // Read 5 parameters
        myfile.read((char*)&a, sizeof(a)); NAUXBASD=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); NCOR=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); NACT=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); NVIR=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); NBF=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);

        // Read MO energy
        EIG= new double[NBF];
        for (int ii=0;ii<NBF;ii++){
            myfile.read((char*)&d, sizeof(d));  EIG[ii]=d;
            myfile.seekg(40-sizeof(d),std::ios::cur);
        }

        // Read B332
        B32size=NAUXBASD*NVIR*NACT;
        B32= new double[B32size];
        for (int iact=0;iact<NACT;iact++){
            for (int ixvrt=0;ixvrt<NAUXBASD*NVIR;ixvrt++){
                myfile.read((char*)&d, sizeof(d));  B32[ixvrt+iact*NAUXBASD*NVIR]=d;
                myfile.seekg(40-sizeof(d),std::ios::cur);
            }
        }

        // Read mp2 corr energy
        myfile.read((char*)&d, sizeof(d));  E2_ref = d;
    }
    myfile.close();
}
