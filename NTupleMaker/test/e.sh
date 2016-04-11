
while read line
do
	#echo "DataBaseDir = DesyTauAnalyses/NTupleMaker/data" >> $line
	#echo "Electron12TriggerEff = DesyTauAnalyses/NTupleMaker/data/LeptonScaleFactors/Electron_Ele12Trigger_eff.root" >> $line
	#echo "Electron12TriggerEff = DesyTauAnalyses/NTupleMaker/data/LeptonScaleFactors/Electron_Ele12Trigger_eff.root" >> $line
	#echo "ElectronIdIsoEff = DesyTauAnalyses/NTupleMaker/data/LeptonScaleFactors/Electron_IdIso0p15_eff.root" >> $line
	#echo "El22LegMC =  hltSingleEle22WP75GsfTrackIsoFilter " >>$line
	#echo "El23LegData =   hltEle23WPLooseGsfTrackIsoFilter  " >>$line
	#echo "SingleMuonTriggerPtCut = 23" >> $line
	#echo "SingleElectronTriggerEtaCut = 2.4">>$line

	#echo "SingleElectronTriggerPtCut = 23" >> $line
	#echo "SingMuonTriggEff = DesyTauAnalyses/NTupleMaker/data/LeptonScaleFactors/Muon_SingleMu_eff.root" >> $line
        # echo	" Muon17TriggerEff = DesyTauAnalyses/NTupleMaker/data/LeptonScaleFactors/Muon_Mu17_eff.root" >>$line
	echo "SingElectronTriggEff = DesyTauAnalyses/NTupleMaker/data/LeptonScaleFactors/Electron_SingleEle_eff.root" >> $line
	#echo "TauFakeRateEff = DesyTauAnalyses/NTupleMaker/data/TFR.root" >> $line

	#echo "MuonIdIsoEff = DesyTauAnalyses/NTupleMaker/data/LeptonScaleFactors/Muon_IdIso0p15_eff.root" >> $line
#	echo "MuonSfDataBarrel = DataMuMu_muonBarrel.root" >> $line 
#	echo "MuonSfDataEndcap = DataMuMu_muonEndcap.root" >> $line
#	echo "MuonSfMcBarrel = DYJetsToMuMu_M-50_muonBarrel.root" >> $line
#	echo "MuonSfMcEndcap = DYJetsToMuMu_M-50_muonEndcap.root" >> $line
done<$1
