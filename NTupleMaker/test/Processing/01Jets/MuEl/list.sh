rm datasets_signal
#Merged_Plots_${signal}400_LSP50_B_OS_np.root

signal=$1

for i in `ls ${signal}*B.root`
do

st=`echo $i   | awk -F "_LSP" '{print $1}' | awk -F "${signal}_" '{print $2}'`
lsp=`echo $i  | awk -F "_LSP" '{print $2}' | awk -F "_" '{print $1}'`

echo $i $st $lsp >> datasets_signal

done

