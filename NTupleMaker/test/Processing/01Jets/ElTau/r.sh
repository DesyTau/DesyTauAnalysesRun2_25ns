

#cp SingleMuon_A.root SingleMuon_InvMuIso_A.root
#cp SingleMuon_B.root SingleMuon_InvMuIso_B.root
#cp SingleMuon_B.root SingleMuon_InvMuIso.root
#cp SingleMuon_C.root SingleMuon_InvMuIso_C.root
#cp SingleMuon_D.root SingleMuon_InvMuIso_D.root

for i `ls *__ext*.root`

do

f=`echo $i | awk -F "__ext[1-4]_" '{print $1}'`
f2=`echo $i | awk -F "_" '{print $2}'`

mv $i ${f}_B_OS.root

done

