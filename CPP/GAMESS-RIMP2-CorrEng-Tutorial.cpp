//  Copyright (C) 2020, Argonne National Laboratory. All Rights Reserved.
//  Licensed under the NCSA open source license

#include "common.h"

double *eij, *eab, *B32;           // They were 2D arrays in Fortran
double *EIG;                       // 1D array
int NAUXBASD,NCOR,NACT,NVIR,NBF;
int B32size;
double E2_ref;

void RIMP2_Energy_Whole_Combined(double *E2);
void Initialization(int argc, char *argv[]);
void run_RIMP2(bool timer_on);
void Finalization();
void Read_Input_File(std::string fname);

int main(int argc, char *argv[]){

    // Read or generate iput data
    Initialization(argc, argv);

    // Run RIMP2
    run_RIMP2(true);

    // Finalize
    Finalization();
    return(0);

}


void run_RIMP2(bool timer_on){

    // Energy
    double E2, E2_diff, Rel_E2_error;

    // Wall time
    double dt;
    double tic,toc;

#if defined(OMP)
    std::cout<<"\n\tRunning the code with OpenMP threading";
#elif defined(OFFLOAD)
    std::cout<<"\n\tRunning the code with OpenMP offloading on GPU";
#else
    std::cout<<"\n\tRunning the code serially";
#endif
#if defined(MKL)
    std::cout<<" with MKL:\n";
#else
    std::cout<<" with a hand-written DGEMM:\n";
#endif


    // Measuing the performance of Correlation Energy Accumulation
    E2 = 0.0;
    tic = timer();
    RIMP2_Energy_Whole_Combined(&E2);
    toc = timer();
    dt = toc-tic;

    // Report the performance data and pass/fail status
    if(timer_on){
      E2_diff = E2 - E2_ref;
      Rel_E2_error = E2_diff/E2_ref;
      Rel_E2_error = (Rel_E2_error > 0) ? (Rel_E2_error) : (-Rel_E2_error);
#if defined(OMP)
      std::cout<<"\t\tNumber of OMP threads                   = "<<omp_get_max_threads()<<"\n";
#endif
//      std::cout<<"\t\tReference MP2 corr. energy              = "<<E2_ref<<"\n";
      std::cout<<"\t\tRel. error of computed MP2 corr. energy = "<<Rel_E2_error<<"\n";
      std::cout<<"\t\tWall time                               = "<<dt<<" sec\n";
      if (Rel_E2_error <= TOL) {
          std::cout<<"\t\tPassed :-) \n";
      }
      else {
          std::cout<<"\t\tFailed :-) \n";
      }
    }

}



void Initialization(int argc, char *argv[]){

    std::string fname;

    // Read input filename
    fname="cor.rand";
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

    // Print out the summary of the input
    std::cout<<"\tNAUXBASD NCOR NACT NVIR NBF = "<<NAUXBASD<<" "<<NCOR<<" "<<NACT<<" "<<NVIR<<" "<<NBF<<"\n";
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


