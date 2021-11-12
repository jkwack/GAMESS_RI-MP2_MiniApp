CPP=icpc
CPP_OFF=icpx

SDIR=./
SRC=${SDIR}/GAMESS-RIMP2-CorrEng-Tutorial.cpp
SRC_V0_CPU=${SDIR}/RIMP2_V0_CPU.cpp
SRC_V1_CPU=${SDIR}/RIMP2_V1_CPU.cpp
SRC_V3_GPU=${SDIR}/RIMP2_V3_GPU.cpp
SRC_V5_GPU=${SDIR}/RIMP2_V5_GPU.cpp

CPPFLAGS_SERIAL=-Ofast -g
CPPFLAGS_SERIAL+= -xCORE-AVX2

CPPFLAGS_OFFLOAD=-fiopenmp -DOFFLOAD -g -O3
CPPFLAGS_OFFLOAD+= -fopenmp-targets=spir64

CPPFLAGS_MKL_OFFLOAD= -fsycl

LDFLAGS_MKL=-qmkl=sequential

LDFLAGS_MKL_OFFLOAD+= -qmkl 

EXE_V0_CPU=rimp2-V0-CPU
EXE_V1_CPU=rimp2-V1-CPU
EXE_V3_GPU=rimp2-V3-GPU
EXE_V5_GPU=rimp2-V5-GPU

EXES=$(EXE_V0_CPU) $(EXE_V1_CPU) $(EXE_V3_GPU) $(EXE_V5_GPU) 

all: $(EXES)

$(EXE_V0_CPU): $(SRC) $(SRC_V0_CPU)
	$(CPP) $(CPPFLAGS_SERIAL) $^ -o $@ 
	rm -rf *.o *.mod *.modmic

$(EXE_V1_CPU): $(SRC) $(SRC_V1_CPU)
	$(CPP) $(CPPFLAGS_SERIAL) -DMKL $^ -o $@ $(LDFLAGS_MKL)
	rm -rf *.o *.mod *.modmic

$(EXE_V3_GPU): $(SRC) $(SRC_V3_GPU)
	$(CPP_OFF) $(CPPFLAGS_OFFLOAD) $(CPPFLAGS_MKL_OFFLOAD) -DMKL $^ -o $@ $(LDFLAGS_MKL_OFFLOAD)
	rm -rf *.o *.mod *.modmic

$(EXE_V5_GPU): $(SRC) $(SRC_V5_GPU)
	$(CPP_OFF) $(CPPFLAGS_OFFLOAD) $(CPPFLAGS_MKL_OFFLOAD) -DMKL $^ -o $@ $(LDFLAGS_MKL_OFFLOAD)
	rm -rf *.o *.mod *.modmic

.SUFFIXES:
.SUFFIXES: .c .o .f90 .cu .cpp .cuf
.PHONY: clean
clean:
	rm -rf *.o *.mod *.modmic $(EXES) $(EXES_MKL)
