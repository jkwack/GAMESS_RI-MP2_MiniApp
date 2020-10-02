FC=ifort
FC_OFF=ifx

SDIR=./
SRC=${SDIR}/GAMESS-RIMP2-CorrEng-Tutorial.f90

EXE_MKL=rimp2-mkl
FFLAGS_MKL=-Ofast -qopenmp -cpp -DCPU -g
LDFLAGS_MKL=-L${MKLROOT}/lib/intel64_lin -I${MKLROOT}/include -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

EXE_MKL_OFFLOAD=rimp2-mkl-offload
FFLAGS_MKL_OFFLOAD=-Ofast -fiopenmp -fopenmp-targets=spir64 -cpp -DINTEL_OFFLOAD -g
LDFLAGS_MKL=-mkl

EXES=$(EXE_MKL) $(EXE_MKL_OFFLOAD)

all: $(EXES)

$(EXE_MKL): $(SRC)
	$(FC) $(FFLAGS_MKL) $^ -o $@ $(LDFLAGS_MKL)
	rm -rf *.o *.mod

$(EXE_MKL_OFFLOAD): $(SRC)
	$(FC_OFF) -mkl $(FFLAGS_MKL_OFFLOAD) $^ -o $@ $(LDFLAGS_MKL)
	rm -rf *.o *.mod


.SUFFIXES:
.SUFFIXES: .c .o .f90 .cu .cpp .cuf
.PHONY: clean
clean:
	rm -rf *.o *.mod *.modmic $(EXES)
