import os 
import numpy as np
import sys

isData=True

cmd = './condorsub_seq_leptau.sh SynchNTupleProducer_tt_2017'


# os.system('bash make_lists.sh')
# os.system('bash make_config.sh')
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

data_samples = ['DATA_TauB','DATA_TauC','DATA_TauD','DATA_TauE','DATA_TauF']

data_conf = 'analysisMacroSynch_tautau_Data.conf'
size=1
if(isData):
    size=len(data_samples)
else:
    size=len(samples)

for i in range(size):
    print('Using '+samples[i])
    if(isData==False):
        Synch(configs[i],samples[i])
    else:
        Synch(data_conf,data_samples[i])

    
