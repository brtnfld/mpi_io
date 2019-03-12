#!/bin/bash -l
#SBATCH -o brtnfld.o%j
#SBATCH -t 00:15:00
#SBATCH --ntasks-per-node=32 # Cori
#SBATCH -C haswell
####SBATCH --ntasks-per-node=24 # Edison
###   1  2  4   8   16  32   64   128
###   32 64 128 256 512 1024 2048 4096
###   1  2  4   8  16  32   64  128  256   512  1024
#### 24 48 96 192 384 768 1536 3072 6144 12288 24576
#SBATCH -N 32
##SBATCH -p regular
#SBATCH -p debug

module load darshan

export MPICH_MPIIO_STATS=1
export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_TIMERS=1
export DARSHAN_DISABLE_SHARED_REDUCTION=1
export DXT_ENABLE_IO_TRACE=4

cd $SCRATCH
WRKDIR=$SCRATCH/$SLURM_JOB_ID
mkdir $WRKDIR
tsk=$SLURM_NTASKS
lfs setstripe -c 12 -S 16m $WRKDIR
cd $WRKDIR
cmdw="a.out"
cleanup=yes
cp $SLURM_SUBMIT_DIR/$cmdw .

#srun -n $tsk ./$cmdw -i
srun -n $tsk ./$cmdw -c

echo $PWD
cd $SLURM_SUBMIT_DIR
if [ -n "$cleanup" ]; then
  rm -r -f $WRKDIR
fi
