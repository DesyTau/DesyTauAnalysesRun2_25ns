import os 
import numpy as np
import sys

cmd = './condorsub_seq_leptau.sh SynchNTupleProducer_tt_2017'

def Synch(config,sample,job='20'):
    os.system(cmd+' '+config+' '+sample+ ' tt' +' '+job)


configs = ['analysisMacroSynch_lept_tt_DY1JetsToLL.conf',
                 'analysisMacroSynch_lept_tt_DY1JetsToLL_ext1.conf',
                 'analysisMacroSynch_lept_tt_DY2JetsToLL.conf',
                 'analysisMacroSynch_lept_tt_DY2JetsToLL_ext1.conf',
                 'analysisMacroSynch_lept_tt_DY3JetsToLL.conf',
                 'analysisMacroSynch_lept_tt_DY3JetsToLL_ext1.conf',
                 'analysisMacroSynch_lept_tt_DY4JetsToLL.conf',
                 'analysisMacroSynch_lept_tt_DYJetsToLL.conf',
                 'analysisMacroSynch_lept_tt_DYJetsToLL_ext1.conf',
                 'analysisMacroSynch_lept_tt_TTTo2L2Nu.conf',
                 'analysisMacroSynch_lept_tt_W1JetsToLNu.conf',
                 'analysisMacroSynch_lept_tt_W2JetsToLNu.conf',
                 'analysisMacroSynch_lept_tt_W3JetsToLNu.conf',
                 'analysisMacroSynch_lept_tt_W4JetsToLNu.conf',
                 'analysisMacroSynch_lept_tt_WJetsToLNu.conf']
samples = ['DY1JetsToLL','DY1JetsToLL_ext1','DY2JetsToLL','DY2JetsToLL_ext1',
                 'DY3JetsToLL','DY3JetsToLL_ext1','DY4JetsToLL','DYJetsToLL',
                 'DYJetsToLL_ext1','TTTo2L2Nu','W1JetsToLNu','W2JetsToLNu',
                 'W3JetsToLNu','W4JetsToLNu','WJetsToLNu']


for i in range(len(configs)):
    print('Using '+samples[i])
    Synch(configs[i],samples[i])

    
