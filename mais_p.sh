#!/bin/sh
#
#This is an example script example.sh
#
#These commands set up the Grid Environment for your job:
#PBS -N opsw1
#PBS -l nice=16,walltime=36:00:00
#PBS -M seuemail@if.uff.br
#PBS -m abe
#PBS -t 0-10

# Gerando sementes
a=$RANDOM
if [ $(($a%2)) -eq 0 ]
then
    a=$(($a+1))
fi

p = ( 0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 )

cd $PBS_O_WORKDIR
./opsw 10000 10 ${p[$PBS_ARRAYID]} $a 0.1 -6 0.1 0.49 100 5000 2500