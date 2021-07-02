#!/bin/csh 
# $1 : nick
# $2 : dataset
# $3 : era
python read_filelist_from_das.py --nick ${1} --query "${2}" --outputfile ${3}/${1}.list
