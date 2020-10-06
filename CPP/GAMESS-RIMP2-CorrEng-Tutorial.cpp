//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include <iostream>
#include <string>
#include <cstdlib>
#include <new>
#include <fstream>

struct rimp2_input {
    double *eij, *eab, *B32;           // They were 2D arrays in Fortran
    double *EIG;                       // 1D array
    int NAUXBASD,NCOR,NACT,NVIR,NBF;
    int NQVV=0;
    int B32size;
    double E2_ref;
} my;

void Initialization(int argc, char *argv[]);


int main(int argc, char *argv[]){

    // Energy
    double E2;

    // Wall time
    double dt;

#if defined(CPU)
    std::cout<<"You are running the code with CPU OpenMP.\n";
#elif defined(INTEL_OFFLOAD)
    std::cout<<"You are running the code with mkl on GPU.\n";
#else
    std::cout<<"You are running the code serially.\n";
#endif

    Initialization(argc, argv);

    std::cout<<"E2_ref ="<<my.E2_ref<<std::endl;
    return(0);

}


void Initialization(int argc, char *argv[]){

    std::string fname;

    // Read input filename
    if (argc > 1) fname=argv[1];
    if (fname=="benz.kern" || fname=="cor.kern" || fname=="c60.kern" || fname=="w30.kern" || fname=="w60.kern") {
        std::cout<<"\tReading data from "<<fname<<std::endl;
        std::cout<<"Read_Input_File needs to be implemented\n";

        my.B32size=my.NAUXBASD*my.NVIR*my.NACT;
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
    std::cout<<"NQVV="<<my.NQVV<<std::endl;

    // Print out the summary of the input
    std::cout<<"\tNAUXBASD NCOR NACT NVIR NBF = "<<my.NAUXBASD<<" "<<my.NCOR<<" "<<my.NACT<<" "<<my.NVIR<<" "<<my.NBF<<"\n";
    std::cout<<"\tNQVV = "<<my.NQVV<<"\n";
    std::cout<<"\tMemory Footprint:\n";
    std::cout<<"\t\tB32[ "<<my.NAUXBASD*my.NVIR<<" , "<<my.NACT<<" ] = "<<my.NAUXBASD*my.NVIR*my.NACT*8.E-6<<" MB\n";
    std::cout<<"\t\teij[ "<<my.NACT<<" , "<<my.NACT<<" ] = "<<my.NACT*my.NACT*8.E-6<<" MB\n";
    std::cout<<"\t\teab[ "<<my.NVIR<<" , "<<my.NVIR<<" ] = "<<my.NVIR*my.NVIR*8.E-6<<" MB\n";
    std::cout<<"\t\tQVV[ "<<my.NVIR<<" , "<<my.NACT<<" , "<<my.NVIR<<" ] = "<<my.NVIR*my.NACT*my.NACT*8.E-6<<" MB\n";

}

