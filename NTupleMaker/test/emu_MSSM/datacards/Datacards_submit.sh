#!/bin/sh 
# $1 - era
# $2 - sample (Data,EMB,DYJets,TTbar,EWK,SMggH,SMothers,bbH,ggH_t,ggH_b,ggH_i,ggA_t,ggA_b,ggA_i)
# $3 - category 
dir=/nfs/dust/cms/user/rasp/Run/emu_MSSM/Feb10

jobname=Datacards_${1}_${2}_${3}
cat > ${dir}/jobs/${jobname}.sh <<EOF1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cd ${CMSSW_BASE}/src
cmsenv
cd -
CreateCards_emu $1 $2 $3 1
EOF1

cat > ${dir}/jobs/${jobname}.submit <<EOF
+RequestRuntime=3000

RequestMemory = 2000

executable = ${dir}/jobs/${jobname}.sh

transfer_executable = True
universe            = vanilla
getenv              = True
Requirements        = OpSysAndVer == "CentOS7"

output              = ${dir}/jobs/${jobname}.out
error               = ${dir}/jobs/${jobname}.error
log                 = ${dir}/jobs/${jobname}.log

queue

EOF
chmod u+x ${dir}/jobs/${jobname}.sh
chmod u+x ${dir}/jobs/${jobname}.submit
condor_submit ${dir}/jobs/${jobname}.submit
