//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include <iostream>
#include <string>
#include <cstdlib>
#include <new>
#include <fstream>
#include <ctime>
#include "mkl.h"

struct rimp2_input {
    double *eij, *eab, *B32;           // They were 2D arrays in Fortran
    double *EIG;                       // 1D array
    int NAUXBASD,NCOR,NACT,NVIR,NBF;
    int NQVV=0;
    int B32size;
    double E2_ref;
} my;

void RIMP2_Energy_Whole_Combined(double *E2);
void Initialization(int argc, char *argv[]);
void Finalization();
void Read_Input_File(std::string fname);



int main(int argc, char *argv[]){

    // Energy
    double E2, E2_diff, Rel_E2_error;

    // Wall time
    double dt;
    std::clock_t tic,toc;

#if defined(CPU)
    std::cout<<"You are running the code with CPU OpenMP.\n";
#elif defined(INTEL_OFFLOAD)
    std::cout<<"You are running the code with mkl on GPU.\n";
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
    tic = std::clock();
    RIMP2_Energy_Whole_Combined(&E2);
    toc = std::clock();
    dt = (toc-tic)/(double)CLOCKS_PER_SEC;


    // Report the performance data and pass/fail status
    E2_diff = E2 - my.E2_ref;
    Rel_E2_error = abs(E2_diff/my.E2_ref);
    std::cout<<"\tResults:\n";
//    std::cout<<"\t\tReference MP2 corr. energy              = "<<my.E2_ref<<"\n";
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
    int nQVV=my.NVIR*my.NQVV*my.NVIR;
    int iQVV;

    QVV = new double[nQVV];

    for(int JACT=0;JACT<my.NACT;JACT++){
        for(int IACTmod=0;IACTmod<=JACT/my.NQVV;IACTmod++){

            // Set a length of blocks in a raw according my.NQVV
            int IACT = IACTmod*my.NQVV;
            if((IACTmod+1)*my.NQVV>JACT+1) { iQVV = JACT - (IACTmod)*my.NQVV + 1; }
            else { iQVV = my.NQVV; }

            // Compute QVV using dgemm
            int m=my.NVIR*iQVV;
            int n=my.NVIR;
            int k=my.NAUXBASD;
            double alpha = 1.0E0;
            double beta = 0.0E0;
            int lda=my.NAUXBASD;
            int ldb=my.NAUXBASD;
            int ldc=my.NVIR*iQVV;
            dgemm("T", "N",
                &m,     &n,     &k,
                &alpha, &my.B32[IACT*my.NAUXBASD*my.NVIR], &lda,
                        &my.B32[JACT*my.NAUXBASD*my.NVIR], &ldb,
                &beta,  QVV,    &ldc);

            // Accumulate E2
            for(int IC=0; IC<iQVV; IC++){
                double E2_t=0.0;
                for(int IB=0; IB<my.NVIR; IB++){
                    for(int IA=0; IA<my.NVIR; IA++){
                        double Tijab = 
                            QVV[IA+ IC*my.NVIR+ IB*my.NVIR*iQVV]
                            / ( my.eij[IACT+IC + JACT*my.NACT]
                                - my.eab[IA + IB*my.NVIR] );
                        double Q_t = 
                            QVV[IA+ IC*my.NVIR+ IB*my.NVIR*iQVV]
                            + QVV[IA+ IC*my.NVIR+ IB*my.NVIR*iQVV];
                        E2_t = E2_t + Tijab * (Q_t - QVV[IB+ IC*my.NVIR+ IA*my.NVIR*iQVV]);
                    }
                }
                double FAC=2.0E0;
                if(IACT+IC == JACT) FAC=1.0E0;
                *E2 = *E2 + FAC*E2_t;
            }
            //std::cout<<IACTmod<<" "<<JACT<<" "<<*E2<<"\n";


        }
    }



    delete[] QVV;
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
            my.NAUXBASD=420;
            my.NCOR=6;
            my.NACT=15;
            my.NVIR=93;
            my.NBF=120;
        }
        else if(fname=="cor.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of cor.kern\n";
            my.NAUXBASD=1512;
            my.NCOR=24;
            my.NACT=54;
            my.NVIR=282;
            my.NBF=384;
        }
        else if (fname=="c60.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of c60.kern\n";
            my.NAUXBASD=3960;
            my.NCOR=60;
            my.NACT=120;
            my.NVIR=360;
            my.NBF=540;
        }
        else if (fname=="w30.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of w30.kern\n";
            my.NAUXBASD=2520;
            my.NCOR=30;
            my.NACT=120;
            my.NVIR=570;
            my.NBF=750;
        }
        else if (fname=="w60.rand") {
            std::cout<<"\tGenerating arbitrary input data with the structure of w60.kern\n";
            my.NAUXBASD=5040;
            my.NCOR=60;
            my.NACT=240;
            my.NVIR=1140;
            my.NBF=1500;
        }
        else{
            std::cout<<"\tError!\n\tOne of the followings should be used as an input:\n";
            std::cout<<"\t\tbenz.kern, cor.kern, c60.kern, w30.kern, or w60.kern for actual data sets, or\n";
            std::cout<<"\t\tbenz.rand, cor.rand, c60.rand, w30.rand, or w60.rand for arbitrary data sets.\n";
        }

        // Generate MO energy
        my.EIG= new double[my.NBF];
        std::fill_n(my.EIG,my.NBF,1);
        for (int ii=0;ii<my.NACT;ii++) {
            my.EIG[ii+my.NCOR] = 2.0;
        }

        // Generate B32
        my.B32size=my.NAUXBASD*my.NVIR*my.NACT;
        my.B32= new double[my.B32size];
        std::fill_n(my.B32,my.B32size,1.0/my.NAUXBASD);

        // Compute the corresponding mp2 corr energy
        my.E2_ref=0.5*(my.NVIR*my.NVIR)*(my.NACT*my.NACT)/(my.NAUXBASD*my.NAUXBASD);
    }

    // Some parameters
    int NOCC = my.NCOR + my.NACT;

    // Virt-Virt MO energy pairs
    my.eab = new double[my.NVIR*my.NVIR];
    for (int IB=0;IB<my.NVIR;IB++){
        for (int IA=0;IA<=IB;IA++){
            my.eab[IA+IB*my.NVIR] = my.EIG[IA+NOCC] + my.EIG[IB+NOCC];
            my.eab[IB+IA*my.NVIR] = my.eab[IA+IB*my.NVIR];
        }
    }

    // occ-occ MO energy pairs
    my.eij = new double[my.NACT*my.NACT];
    for (int JJ=0;JJ<my.NACT;JJ++){
        for (int II=0;II<=JJ;II++){
            my.eij[II+JJ*my.NACT] = my.EIG[II+my.NCOR] + my.EIG[JJ+my.NCOR];
        }
    }

    // Read the second command line argument for NQVV
    if (argc > 2) {my.NQVV=std::atoi(argv[2]);}
    else {my.NQVV=my.NACT;}
//    std::cout<<"NQVV="<<my.NQVV<<std::endl;

    // Print out the summary of the input
    std::cout<<"\tNAUXBASD NCOR NACT NVIR NBF = "<<my.NAUXBASD<<" "<<my.NCOR<<" "<<my.NACT<<" "<<my.NVIR<<" "<<my.NBF<<"\n";
    std::cout<<"\tNQVV = "<<my.NQVV<<"\n";
    std::cout<<"\tMemory Footprint:\n";
    std::cout<<"\t\tB32[ "<<my.NAUXBASD*my.NVIR<<" , "<<my.NACT<<" ] = "<<my.NAUXBASD*my.NVIR*my.NACT*8.E-6<<" MB\n";
    std::cout<<"\t\teij[ "<<my.NACT<<" , "<<my.NACT<<" ] = "<<my.NACT*my.NACT*8.E-6<<" MB\n";
    std::cout<<"\t\teab[ "<<my.NVIR<<" , "<<my.NVIR<<" ] = "<<my.NVIR*my.NVIR*8.E-6<<" MB\n";
    std::cout<<"\t\tQVV[ "<<my.NVIR<<" , "<<my.NACT<<" , "<<my.NVIR<<" ] = "<<my.NVIR*my.NACT*my.NVIR*8.E-6<<" MB\n";

}



void Finalization(){
    delete[] my.B32;
    delete[] my.EIG;
    delete[] my.eij;
    delete[] my.eab;
}



void Read_Input_File(std::string fname){
    std::ifstream myfile;
    int a;
    double d;

    myfile.open(fname,std::ios::in | std::ios::binary);

    if (myfile.is_open()){

        // Read 5 parameters
        myfile.read((char*)&a, sizeof(a)); my.NAUXBASD=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); my.NCOR=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); my.NACT=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); my.NVIR=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);
        myfile.read((char*)&a, sizeof(a)); my.NBF=a;
        myfile.seekg(40-sizeof(a),std::ios::cur);

        // Read MO energy
        my.EIG= new double[my.NBF];
        for (int ii=0;ii<my.NBF;ii++){
            myfile.read((char*)&d, sizeof(d));  my.EIG[ii]=d;
            myfile.seekg(40-sizeof(d),std::ios::cur);
        }

        // Read B332
        my.B32size=my.NAUXBASD*my.NVIR*my.NACT;
        my.B32= new double[my.B32size];
        for (int iact=0;iact<my.NACT;iact++){
            for (int ixvrt=0;ixvrt<my.NAUXBASD*my.NVIR;ixvrt++){
                myfile.read((char*)&d, sizeof(d));  my.B32[ixvrt+iact*my.NAUXBASD*my.NVIR]=d;
                myfile.seekg(40-sizeof(d),std::ios::cur);
            }
        }

        // Read mp2 corr energy
        myfile.read((char*)&d, sizeof(d));  my.E2_ref = d;
    }
    myfile.close();
}
