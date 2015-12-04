#!/bin/csh
./hadd.csh SingleMuon_2015D_05Oct
./hadd.csh SingleMuon_2015D_PRv4
rm SingleMuon_2015D.root
hadd SingleMuon_2015D.root SingleMuon_2015D_05Oct.root SingleMuon_2015D_PRv4.root
