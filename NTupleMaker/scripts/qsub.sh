#!/bin/sh 
# $1 - code name
# $2 - config file
# $3 - sample name
# $4 - channel [optional]
cat > $3.zsh <<EOF
#!/bin/zsh
#
#(make sure the right shell will be used)
#$ -S /bin/zsh
#
#(the cpu time for this job)
#$ -l h_rt=23:59:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=2G
#
#(use hh site)
#$ -l site=hh 
#(stderr and stdout are merged together to stdout)
#$ -j y
#
# use SL6
#$ -l os=sld6
#
# use current dir and current environment
#$ -cwd
#$ -V
#
#$ -o $3.out
#
#$ -e $3.err
$@

EOF

chmod u+x $3.zsh
qsub $3.zsh
