#!/usr/bin/env python

import commands
import re
import os

import sys
sys.path.append('./')

###########################################
###########################################

def preAnalyze( ana, sample, runInSeries=False):

    #os.system( 'mkdir batch/' )

    stream = "TauTau"
    print "Stream ", stream

    #os.system('python makeTreeSkimmerTauTau_Summer14.py')

    if runInSeries:
         print "Running in series via the command ..."
         configFileName = "./Configs/preAnalyzerTauTau_Summer15_%s_%s_cfg.py" % (sample,ana)
         runInSeries   = 'preAnalyzer'+stream+' '+configFileName+'\n'
         print runInSeries
         os.system(runInSeries)

    else:
        nameJob = 'job_'+sample+'_'+ana+'_'+stream                                                                                                                           
        if 'Data' in sample:
            nameJob = 'job_'+sample+'_'+stream
        else:
            nameJob = 'job_'+sample+'_'+stream+'_'+ana
        fileJob = 'batch/'+nameJob+'.sh'
        fileLog = 'batch/log/'+nameJob+'.txt'
        print "Creating the shell file : "+nameJob
        
        ##
        f = open(fileJob,'w')
        f.write('#!/bin/sh\n\n')
        configFileName = "./Configs/preAnalyzerTauTau_Summer15_%s_%s_cfg.py" % (sample,ana)
        f.write('cfg='+configFileName+'\n')
        f.write('cat > $cfg.zsh <<EOF\n')
        f.write('#!/bin/zsh\n')
        f.write('#\n')
        f.write('#(make sure the right shell will be used)\n')
        f.write('#$ -S /bin/zsh\n')
        f.write('#\n')
        f.write('#(the cpu time for this job)\n')
        f.write('#$ -l h_cpu=02:29:00\n')
        f.write('#\n')
        f.write('#(the maximum memory usage of this job)\n')
        f.write('#$ -l h_vmem=2000M\n')
        f.write('#\n')
        f.write('#(use hh site)\n')
        f.write('#$ -l site=hh\n')
        f.write('#(stderr and stdout are merged together to stdout)\n')
        f.write('#$ -j y\n')
        f.write('#\n')
        f.write('# use SL6\n')
        f.write('#$ -l os=sld6\n')
        f.write('#\n')
        f.write('# use current dir and current environment\n')
        f.write('#$ -cwd\n')
        f.write('#$ -V\n')
        f.write('#\n')
        f.write('#$ -o $cfg.out\n')
        f.write('#\n')
        f.write('#$ -e $cfg.err\n')
        f.write('preAnalyzer'+stream+' $cfg\n')
        
        f.write('EOF\n')

        f.write('rm $cfg.out\n')
        f.write('rm $cfg.err\n')
        f.write('chmod u+x $cfg.zsh\n')
        f.write('qsub $cfg.zsh\n')
        f.close()
        os.system('chmod u+x batch/*.sh')

os.system('python preAnalyzerTauTau_Summer15.py')
    
###########################################
###########################################
#Data
preAnalyze("nominal","Tau_Run2015D_05Oct2015",False)
preAnalyze("nominal","Tau_Run2015D_PRv4",False)

analysis = {"nominal", "TauUp", "TauDown"}
for ana in analysis:
#Signal
    preAnalyze(ana,"GGFH125",False)
#preAnalyze(ana,"VBFH125",False)
    preAnalyze(ana,"SUSYGGH160",False)
#preAnalyze(ana,"SUSYBBH160",False)
    
# Background
    preAnalyze(ana,"DYJets_TauTau",False)
    preAnalyze(ana,"WJetsToLNu",False)
    preAnalyze(ana,"TTJets",False)
    preAnalyze(ana,"SingleTop_t",False)
    preAnalyze(ana,"SingleAntiTop_t",False)
    preAnalyze(ana,"SingleTop_tW",False)
    preAnalyze(ana,"SingleAntiTop_tW",False)
    preAnalyze(ana,"VVTo2L2Nu",False)
    preAnalyze(ana,"WZJets",False)
    preAnalyze(ana,"WWToLNuQQ",False)
    preAnalyze(ana,"WZTo1L1Nu2Q",False)
    preAnalyze(ana,"WZTo1L3Nu",False)
    preAnalyze(ana,"WZTo2L2Q",False)
    preAnalyze(ana,"ZZTo2L2Q",False)
    preAnalyze(ana,"ZZTo4L",False)
