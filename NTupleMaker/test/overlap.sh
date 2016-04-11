
while read line
	do

		RunCor sm ${line}_B_OS_np Ranking
	done<$1
