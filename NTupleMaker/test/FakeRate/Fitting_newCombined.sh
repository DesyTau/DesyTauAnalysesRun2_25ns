#!/bin/sh
#cd /nfs/dust/cms/user/ywen/CMSSW/combineTool_New/CMSSW_8_1_0/src/
cd /nfs/dust/cms/user/cardinia/Combine/CMSSW_8_1_0/src/
eval `scramv1 runtime -sh`
cd -

templateFittingETauFR ETauFRAgainstEleVLooseLt1p460.root | tee ./coutFR/ETauFRCoutVLooseLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.03889" ETauFRAgainstEleVLooseLt1p460.txt -o  WorkSpaceVLooseLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVLooseLt1p460.root | tee ./coutFR/ScaleVLooseLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VLooseLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRVLooseLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVLooseLt1p460.txt -w WorkSpaceVLooseLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVLooseLt1p460.root -m 90 -o impacts_VLooseLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VLooseLt1p460.json -o ./coutFR/impacts_VLooseLt1p460


templateFittingETauFR ETauFRAgainstEleLooseLt1p460.root | tee ./coutFR/ETauFRCoutLooseLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00752276" ETauFRAgainstEleLooseLt1p460.txt -o WorkSpaceLooseLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceLooseLt1p460.root | tee ./coutFR/ScaleLooseLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/LooseLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRLooseLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleLooseLt1p460.txt -w WorkSpaceLooseLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceLooseLt1p460.root -m 90 -o impacts_LooseLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_LooseLt1p460.json -o ./coutFR/impacts_LooseLt1p460


templateFittingETauFR ETauFRAgainstEleMediumLt1p460.root | tee ./coutFR/ETauFRCoutMediumLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00196609" ETauFRAgainstEleMediumLt1p460.txt -o WorkSpaceMediumLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceMediumLt1p460.root | tee ./coutFR/ScaleMediumLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/MediumLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRMediumLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleMediumLt1p460.txt -w WorkSpaceMediumLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceMediumLt1p460.root -m 90 -o impacts_MediumLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_MediumLt1p460.json -o ./coutFR/impacts_MediumLt1p460


templateFittingETauFR ETauFRAgainstEleTightLt1p460.root | tee ./coutFR/ETauFRCoutTightLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.000822007" ETauFRAgainstEleTightLt1p460.txt -o WorkSpaceTightLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceTightLt1p460.root | tee ./coutFR/ScaleTightLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/TightLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRTightLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleTightLt1p460.txt -w WorkSpaceTightLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceTightLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceTightLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceTightLt1p460.root -m 90 -o impacts_TightLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_TightLt1p460.json -o ./coutFR/impacts_TightLt1p460


templateFittingETauFR ETauFRAgainstEleVTightLt1p460.root | tee ./coutFR/ETauFRCoutVTightLt1p460.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.000429945" ETauFRAgainstEleVTightLt1p460.txt -o WorkSpaceVTightLt1p460.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVTightLt1p460.root | tee ./coutFR/ScaleVTightLt1p460.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VTightLt1p460Pull.txt
PostFitShapesFromWorkspace -o ETauFRVTightLt1p460_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVTightLt1p460.txt -w WorkSpaceVTightLt1p460.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightLt1p460.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightLt1p460.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVTightLt1p460.root -m 90 -o impacts_VTightLt1p460.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VTightLt1p460.json -o ./coutFR/impacts_VTightLt1p460


templateFittingETauFR ETauFRAgainstEleVLooseGt1p558.root | tee ./coutFR/ETauFRCoutVLooseGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.0991319" ETauFRAgainstEleVLooseGt1p558.txt -o  WorkSpaceVLooseGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVLooseGt1p558.root | tee ./coutFR/ScaleVLooseGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VLooseGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRVLooseGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVLooseGt1p558.txt -w WorkSpaceVLooseGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVLooseGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVLooseGt1p558.root -m 90 -o impacts_VLooseGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VLooseGt1p558.json -o ./coutFR/impacts_VLooseGt1p558


templateFittingETauFR ETauFRAgainstEleLooseGt1p558.root | tee ./coutFR/ETauFRCoutLooseGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.0192047" ETauFRAgainstEleLooseGt1p558.txt -o WorkSpaceLooseGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVERq_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceLooseGt1p558.root | tee ./coutFR/ScaleLooseGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/LooseGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRLooseGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleLooseGt1p558.txt -w WorkSpaceLooseGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceLooseGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceLooseGt1p558.root -m 90 -o impacts_LooseGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_LooseGt1p558.json -o ./coutFR/impacts_LooseGt1p558


templateFittingETauFR ETauFRAgainstEleMediumGt1p558.root | tee ./coutFR/ETauFRCoutMediumGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00366973" ETauFRAgainstEleMediumGt1p558.txt -o WorkSpaceMediumGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceMediumGt1p558.root | tee ./coutFR/ScaleMediumGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/MediumGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRMediumGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleMediumGt1p558.txt -w WorkSpaceMediumGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceMediumGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceMediumGt1p558.root -m 90 -o impacts_MediumGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_MediumGt1p558.json -o ./coutFR/impacts_MediumGt1p558


templateFittingETauFR ETauFRAgainstEleTightGt1p558.root | tee ./coutFR/ETauFRCoutTightGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.00126043" ETauFRAgainstEleTightGt1p558.txt -o WorkSpaceTightGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceTightGt1p558.root | tee ./coutFR/ScaleTightGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/TightGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRTightGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleTightGt1p558.txt -w WorkSpaceTightGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceTightGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceTightGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceTightGt1p558.root -m 90 -o impacts_TightGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_TightGt1p558.json -o ./coutFR/impacts_TightGt1p558


templateFittingETauFR ETauFRAgainstEleVTightGt1p558.root | tee ./coutFR/ETauFRCoutVTightGt1p558.txt
text2workspace.py -m 90 -P HiggsAnalysis.KITHiggsToTauTau.datacards.zttmodels:ztt_eff --PO "eff=0.000626589" ETauFRAgainstEleVTightGt1p558.txt -o WorkSpaceVTightGt1p558.root
combine -m 90  -M FitDiagnostics --robustFit=1 --preFitValue=1. --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1" -n "" WorkSpaceVTightGt1p558.root | tee ./coutFR/ScaleVTightGt1p558.txt
./compare.py -a fitDiagnostics.root | tee ./coutFR/VTightGt1p558Pull.txt
PostFitShapesFromWorkspace -o ETauFRVTightGt1p558_PostFitShape.root -m 90 -f fitDiagnostics.root:fit_s --postfit --sampling --print -d ETauFRAgainstEleVTightGt1p558.txt -w WorkSpaceVTightGt1p558.root
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightGt1p558.root --doInitialFit --robustFit 1 --redefineSignalPOIs r --rMin=0.3 --rMax=2.5 --cminFallbackAlgo "Minuit2,0:1"
combineTool.py -M Impacts -m 90 -d WorkSpaceVTightGt1p558.root --doFits --robustFit 1 --redefineSignalPOIs r --cminFallbackAlgo "Minuit2,0:1" --parallel 5
combineTool.py -M Impacts -d WorkSpaceVTightGt1p558.root -m 90 -o impacts_VTightGt1p558.json --redefineSignalPOIs r
plotImpacts.py -i impacts_VTightGt1p558.json -o ./coutFR/impacts_VTightGt1p558

