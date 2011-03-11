#!/bin/sh

if [ $# -ne 11 ]
then
    echo 'Use '$0 '<n> <K> <p-SW> <semente> <sem_rede> <eps> <J> <q> <mcs> <intervalo para salvar>'
    echo 'As saídas são arquivos .dat para a série temporal, .bin para o estado de número aleatório e .st para o estado inicial.'
    exit
fi

echo 'Compilando...'
icc zero2_gsl.c -o zero.out -O3 -ipo -xT -parallel -lgsl -lgslcblas -lm
icc opinionsw.c -o opsw.out -O3 -ipo -xT -parallel -lgsl -lgslcblas -lm

# Para construtir estado
n=$1
K=$2
p=$3
sem=$4
sem_rede=$5
eps=$6
J=$7
q=$8
rho=$9

# Para executar simulação
Tmcs=${10}
Dmcs=${11}

# formando o nome de saida
saida='N='$n'_p='$p'_rho_0='$rho
estado=$saida.st
binsem=$saida.bin
ts=$saida.dat

echo 'Executando estado inicial'
./zero.out $n $K $p $sem $sem_rede $eps $J $q $estado $binsem $rho
echo 'Executando dinamica de opiniao'
./opsw.out $estado $Tmcs $Dmcs $ts $binsem

