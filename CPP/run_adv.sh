#!/bin/bash

# This is an example of how to run advisor. As input, it takes in what
# you want to run.
#
# Example:
# $ ./run_adv.sh ./rimp2-CPP-V45-MKL-OFFLOAD-OPT cor.rand

if [ $# -eq 0 ]; then
    echo "To run, you need to pass executable and arguments you want, like:"
    echo ""
    echo "$ ./run_adv.sh ./rimp2-CPP-V45-MKL-OFFLOAD-OPT cor.rand"
    exit 1
fi

set -x

advixe-cl --collect=roofline --profile-gpu --project-dir=./out -data-limit=0 --search-dir src:r=./ -- $@

advixe-cl --report=roofline --gpu --project-dir=./out --report-output=./roofline.out.html

advixe-python ${ADVISOR_2021_DIR}/pythonapi/examples/survey_gpu.py ./out > ./advisor_raw_data_out.txt
