#rm corrupted_jobs

qstat > jobs
cat jobs | awk -F " " '{print $1}' > j

while read line 
do

echo cat *$line* | head -3 | tail -1 >> corrupted_jobs
done<j

