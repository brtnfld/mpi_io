#!/bin/bash -l
#SBATCH -o brtnfld.o%j
#SBATCH -t 00:20:00
#SBATCH --ntasks-per-node=32 # Cori
#SBATCH -C haswell 
##SBATCH --ntasks-per-node=24 # Edison
#SBATCH -n 16384 #8192 #4096
#SBATCH -p regular
##SBATCH -p debug
cd $SCRATCH
WRKDIR=$SCRATCH/$SLURM_JOB_ID
mkdir $WRKDIR
tsk=$SLURM_NTASKS #16384 #8192 # 4096
lfs setstripe -c 12 -S 16m $WRKDIR
cd $WRKDIR
cmd="a.out"
#script="bench.sh"
##cleanup=yes
cp $SLURM_SUBMIT_DIR/$cmd .
#cp $SLURM_SUBMIT_DIR/$script .
for i in `seq 1 10`; do
  srun -n $tsk ./$cmd
done
for i in `seq 1 10`; do
  srun -n $tsk ./$cmd 0
done
#cp sum_* $SLURM_SUBMIT_DIR/
#cp file* $SLURM_SUBMIT_DIR/
cd $SLURM_SUBMIT_DIR
if [ -n "$cleanup" ]; then
  rm -r -f $WRKDIR
fi
