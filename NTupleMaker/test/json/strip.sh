
sed 's/]],/&\n/g' d.json  > temp; sed -e 's/{//g' -e 's/}//g' -e 's/"//g'  temp > t ; mv t temp
sed -i 's/^[ \t]*//' temp
cat temp


