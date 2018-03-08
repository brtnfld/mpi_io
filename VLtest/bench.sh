#!/bin/bash

HOST=$(uname -n)
RUN=''
if [ "$HOST" == 'jelly' ];then
    RUN="$HOME/packages/nocache/nocache"
fi

nt=10
##for j in 4096 8192 16384 32768 65536 131072 262144
for j in 4096 8192 16384 32768 65536 131072 262144
do
 rm -f VL_timing.txt
 for i in `seq 1 $nt`
  do
    ## echo $RUN a.out $j
    $RUN a.out $j
    sleep 5
  done
  perl -lane 'for $c (0..$#F){$t[$c] += $F[$c]}; END{for $c (0..$#t){print $t[$c]/$.}}' VL_timing.txt >> $j
done
