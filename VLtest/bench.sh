#!/bin/bash

HOST=$(uname -n)
RUN=''
if [ "$HOST" == 'jelly' ];then
    RUN="$HOME/packages/nocache/nocache"
fi
cmd=$1
nt=10
##for j in 4096 8192 16384 32768 65536 131072 262144
for j in 4096 8192 16384 32768 65536 131072 262144
do
 rm -f VL_timing.txt
 for i in `seq 1 $nt`
  do
    $RUN ./$cmd -w $j
    if [ "$i" -eq "1" ] ; then wc -c < h5ex_t_vlen.h5 > filesize_$j ; fi
    $RUN ./$cmd -r $j
#    $RUN ./$cmd -wv $j
#    if [ "$i" -eq "1" ] ; then wc -c < h5ex_t_vlen.h5 >> filesize_$j ; fi
#    $RUN ./$cmd -rv $j
  done
  perl -lane 'for $c (0..$#F){$t[$c] += $F[$c]}; END{for $c (0..$#t){print $t[$c]/$.}}' VL_timing.txt >> sum_$j
done
