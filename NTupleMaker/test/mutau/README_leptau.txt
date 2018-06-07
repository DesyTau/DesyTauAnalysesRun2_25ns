Instructions on running the mutau channel (the other channels can be
run in a similar way).

1. use the configuration files $analysisMacroSynch_lept_mt_MC17.conf
and $analysisMacroSynch_lept_mt_DATA_SingleMuon.conf for the MC samples
and the data requiring the single muon trigger (the Xtrigger is
implemented as well in this config file).

2. run $make_config.sh to create the configuration files for each MC
samples. The macro assigns to each sample the corresponding pileup histogram.

3. edit and run $make_lists.sh to create the file lists for the MC
samples and data

4. use the $run_mt.sh script to run the SynchNTupleProducer_2017 macro
on all the MC samples.

This script requires having the $condorsub_seq_leptau.sh and the
$HTC_submit.sh macros in the same area where the code is run.

Note: the uploaded make_config.sh macro is intended for the mutau
channel, just change "_mt" with "_et" or "_tt" to run on the etau or
tau tau channels. 
