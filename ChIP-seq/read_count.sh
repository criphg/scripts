for a in `ls Q30*_sorted.bed`; do
	tag=`echo $a | awk -F"_" '{print$2}'`;
	echo "Start: $a ,tag = $tag";
	(for i in `ls *_region*`; do
		echo $i; 
		chrm=`cat $i | awk '{print$1}'`;
		echo "/"$chrm/"";
		awk -v var=$chrm '{if($1 == var) print$0}' $a > Chr_temp.bed; 
		(echo $i;bedtools intersect -a $i -b Chr_temp.bed|wc -l) >> $tag"_CountRead.csv" ; 
	done;)
#	rm -f Chr_temp.bed;
done;
