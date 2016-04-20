#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#
#(the cpu time for this job)
#$ -l h_cpu=1:59:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=1000M
#
#(use hh site)
#$ -l site=hh
#(stderr and stdout are merged together to stdout)
#$ -j y
#
# use SL5
#$ -l os=sld6
#
# use current dir and current environment
#$ -cwd
#$ -V


channel=$2 
era=$3

era=25ns
#era=InvMuIso
#era=eltau

channel=mutau
#channel=eltau


while read line
do


if [ ! -f $era/${line}_B_OS.root ] 
then
SUSY$channel analysisMacroSUSY_MC_B.conf ${line} $era

fi

done<$1



