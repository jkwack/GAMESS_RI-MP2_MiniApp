CPP=icpc
CPP_OFF=icpx

SDIR=./
SRC=${SDIR}/GAMESS-RIMP2-CorrEng-Tutorial.cpp

EXE_MKL=rimp2-CPP
FFLAGS_MKL=-Ofast -g
LDFLAGS_MKL=-mkl=sequential

EXE_MKL_OMP=rimp2-CPP-OMP
FFLAGS_MKL_OMP=-Ofast -qopenmp -DOMP -g
LDFLAGS_MKL_OMP=-mkl=sequential

EXE_MKL_OFFLOAD=rimp2-CPP-OFFLOAD
FFLAGS_MKL_OFFLOAD=-Ofast -fiopenmp -fopenmp-targets=spir64 -DOFFLOAD -g
LDFLAGS_MKL_OFFLOAD=-mkl

EXES=$(EXE_MKL) $(EXE_MKL_OMP) $(EXE_MKL_OFFLOAD)

all: $(EXES)

$(EXE_MKL): $(SRC)
	$(CPP) $(FFLAGS_MKL) $^ -o $@ $(LDFLAGS_MKL)
	rm -rf *.o *.mod *.modmic

$(EXE_MKL_OMP): $(SRC)
	$(CPP) $(FFLAGS_MKL_OMP) $^ -o $@ $(LDFLAGS_MKL_OMP)
	rm -rf *.o *.mod *.modmic

$(EXE_MKL_OFFLOAD): $(SRC)
	$(CPP_OFF) $(FFLAGS_MKL_OFFLOAD) $^ -o $@ $(LDFLAGS_MKL_OFFLOAD)
	rm -rf *.o *.mod *.modmic


.SUFFIXES:
.SUFFIXES: .c .o .f90 .cu .cpp .cuf
.PHONY: clean
clean:
	rm -rf *.o *.mod *.modmic $(EXES)
