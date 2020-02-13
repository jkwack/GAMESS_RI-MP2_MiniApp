#!/bin/bash
make clean
make all
NMPI=2 INPUT=$1 NQVV=$2 ./run_gpu.sh
NMPI=2 INPUT=$1 NQVV=$2 ./run_cpu.sh


