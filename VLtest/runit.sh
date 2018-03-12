#!/bin/bash -l
##SBATCH -p regular
#SBATCH -p debug
#SBATCH -t 00:30:00
#SBATCH -N 1 #171 #1366 #1366 #171 #342
#SBATCH -o brtnfld.o%j
##SBATCH -C haswell   #Use Haswell nodes
#Edison has 24 cores per node
cd $SCRATCH
WRKDIR=$SCRATCH/$SLURM_JOB_ID
mkdir $WRKDIR
lfs setstripe -c 12 -S 16m $WRKDIR
cd $WRKDIR
cmd="a.out"
script="bench.sh"
cleanup=yes
cp $SLURM_SUBMIT_DIR/$cmd .
cp $SLURM_SUBMIT_DIR/$script .
./$script
cp sum_* $SLURM_SUBMIT_DIR/
cd $SLURM_SUBMIT_DIR
if [ -n "$cleanup" ]; then
  rm -r -f $WRKDIR
fi
