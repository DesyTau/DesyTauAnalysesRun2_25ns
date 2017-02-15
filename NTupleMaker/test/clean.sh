for i in `ls job*`
do
rm $i
done
rm *.o*
rm *sh.e*
rm -fr dir_*
rm list_*
rm pmu*
rm pel*
rm pWte*
rm *~
for i in `ls Jobs/job*`
do
	rm Jobs/$i*
done
