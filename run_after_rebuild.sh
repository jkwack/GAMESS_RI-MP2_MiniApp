#!/bin/bash
if [ "$1" != "" ]; then
  input=$1
else
  input=benz.kern
fi
make clean
make all
./run_gpu.sh $input $2
./run_cpu.sh $input $2


