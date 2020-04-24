#!/bin/sh 
# $1 - executable
# $2 - filelist 

cat > $2.submit <<EOF
+RequestRuntime=10000

RequestMemory = 2000

executable = $1

transfer_executable = True
universe            = vanilla
getenv              = True
Requirements        = OpSysAndVer == "SL6"

output              = $2.out
error               = $2.error
log                 = $2.log

queue

EOF

chmod u+x $2.submit
condor_submit $2.submit