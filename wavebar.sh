#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -c 1# specify 6 threads per process
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err # optional, name of the stderr, using job and first node values
#SBATCH -J beji_ale_rbim_rk4

date

module purge
module load intel/19.0.5

export HOME_DIR=/home/$USER
export WORK_DIR=/work/$USER

##SLURM_JOBID: Job ID number given to this job
##SLURM_JOB_NODELIST: List of nodes allocated to the job
##SLURM_SUBMIT_DIR: Directory where the sbatch command was executed

mkdir -p $WORK_DIR/$SLURM_JOB_NAME.$SLURM_JOBID
cd $WORK_DIR/$SLURM_JOB_NAME.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.sh .
cp $SLURM_SUBMIT_DIR/*.ipt .
cp $SLURM_SUBMIT_DIR/*.f90 .
cp $SLURM_SUBMIT_DIR/*.out .
cp $SLURM_SUBMIT_DIR/*.err .
ifort -qopenmp wavebar.f90 -mkl

./a.out

# Mark the time it finishes.
date
# exit the job
exit 0
