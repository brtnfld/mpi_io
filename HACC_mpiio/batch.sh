#!/bin/bash
#BSUB -P CSC298
#BSUB -W 00:25
# power 42
#BSUB -nnodes 64 
#BSUB -J mpiio 
#BSUB -o mpiio.%J
#BSUB -e mpiio.%J
##SMT1 -- 1 HW Thread per physical core
##SMT4 -- All 4 HW threads are active (Default)
##BSUB -alloc_flags smt1
# 42 physical cores, (21 each cpu), per node
# 84 per cpu, 168 total

##BSUB -alloc_flags maximizegpfs

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export MPICH_MPIIO_STATS=1
export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_TIMERS=1
export DARSHAN_DISABLE_SHARED_REDUCTION=1
export DXT_ENABLE_IO_TRACE=1

EXEC=Sample_mpio_measure_time

JID=$LSB_JOBID
cd $MEMBERWORK/csc298
mkdir mpiio.$JID
cd mpiio.$JID
cp $LS_SUBCWD/$EXEC .

NPROCS="42 84 168 336 672 1344 2688 "
for i in ${NPROCS}
do
  for j in {1..10}
  do
    jsrun -n $i ./$EXEC -i
  done 
  ls -aolF
done
