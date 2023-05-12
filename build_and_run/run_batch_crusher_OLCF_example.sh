#!/bin/bash
#SBATCH -A CHM135_crusher
#SBATCH -J GAMESS_RIMP2
#SBATCH -o %x-%j.out
#SBATCH -t 02:00:00
#SBATCH --threads-per-core=2
#SBATCH --core-spec=0
#SBATCH -N 16

cd /gpfs/alpine/scratch/tiwari/chm135/apps/mini-apps/GAMESS_RI-MP2_MiniApp

source source_me_crusher_OLCF

inputs='w30.rand w60.rand cor.rand c60.rand'

NNODES=$SLURM_JOB_NUM_NODES

if [ "$NNODES" == "1" ]; then
	# 1 node
	iMPIs='1 2 4 8'
elif [ "$NNODES" == "2" ]; then
	# 2 nodes
	iMPIs='16'
elif [ "$NNODES" == "4" ]; then
	# 4 nodes
	iMPIs='32'
elif [ "$NNODES" == "8" ]; then
	# 8 nodes
	iMPIs='64'
elif [ "$NNODES" == "16" ]; then
	#16 nodes
	iMPIs='128'
fi

run_gpu=0
if [ "$run_gpu" == "1" ]; then
for iINPUT in ${inputs}; 
do
  for iMPI in ${iMPIs}; 
  do
    output_file="${SLURM_JOBID}.GPU.NNODES=${NNODES}.NMPI=${iMPI}.INPUT=${iINPUT}"
    set -x
    NMPI=$iMPI INPUT=$iINPUT ./run_gpu.sh > ${output_file} 2>&1
    # post-process the output
    maximum=`cat ${output_file} | grep "Wall time (maximum)" | awk '{print $5}'`
    minimum=`cat ${output_file} | grep "Wall time (minimum)" | awk '{print $5}'`
    mean=`cat ${output_file} | grep "Wall time (mean)" | awk '{print $5}'`
    walltime=`cat ${output_file} | grep ^real | awk '{print $2}'`
    echo ${NNODES} ${SLURM_JOBID} ${iINPUT} ${iMPI} $minimum $mean $maximum $walltime >> global_GPU_NNODES=${NNODES}_${SLURM_JOBID}.log
    set +x
  done
done
fi

run_cpu=1
if [ "$run_cpu" == "1" ]; then
for iINPUT in ${inputs};
do
  iMPI=$((NNODES*16))
  output_file="${SLURM_JOBID}.CPU.NNODES=${NNODES}.NMPI=${iMPI}.INPUT=${iINPUT}"
  set -x
  # run the job_stepper

  NMPI=$iMPI INPUT=$iINPUT ./run_cpu.sh > ${output_file} 2>&1
  set +x

  # post-process the output
  maximum=`cat ${output_file} | grep "Wall time (maximum)" | awk '{print $5}'`
  minimum=`cat ${output_file} | grep "Wall time (minimum)" | awk '{print $5}'`
  mean=`cat ${output_file} | grep "Wall time (mean)" | awk '{print $5}'`
  walltime=`cat ${output_file} | grep ^real | awk '{print $2}'`
  echo ${NNODES} ${SLURM_JOBID} ${iINPUT} ${iMPI} $minimum $mean $maximum $walltime >> global_CPU_NNODES=${NNODES}_${SLURM_JOBID}.log
done
fi




