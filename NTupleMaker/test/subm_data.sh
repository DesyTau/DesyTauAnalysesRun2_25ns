

cat ${2}/${1} | grep "2016B" > ${2}/${1}B
echo ${1}B > ${1}B
qsub -l h_vmem=1500M -l h_rt=0:45:00 run_data.sh ${1}B ${2} grid


cat ${2}/${1} | grep "2016C" > ${2}/${1}C
echo ${1}C > ${1}C
qsub -l h_vmem=1500M -l h_rt=0:45:00 run_data.sh ${1}C ${2} grid

cat ${2}/${1} | grep "2016D" > ${2}/${1}D
echo ${1}D > ${1}D
qsub -l h_vmem=1500M -l h_rt=0:45:00 run_data.sh ${1}D ${2} grid

cat ${2}/${1} | grep "2016E" > ${2}/${1}E
echo ${1}E > ${1}E
qsub -l h_vmem=1500M -l h_rt=0:45:00 run_data.sh ${1}E ${2} grid


cat ${2}/${1} | grep "2016F" > ${2}/${1}F
echo ${1}F > ${1}F
qsub -l h_vmem=1500M -l h_rt=0:45:00 run_data.sh ${1}F ${2} grid

cat ${2}/${1} | grep "2016G" > ${2}/${1}G
echo ${1}G > ${1}G
qsub -l h_vmem=1500M -l h_rt=0:45:00 run_data.sh ${1}G ${2} grid

cat ${2}/${1} | grep "2016H" > ${2}/${1}H
echo ${1}H > ${1}H
qsub -l h_vmem=1500M -l h_rt=0:45:00 run_data.sh ${1}H ${2} grid

