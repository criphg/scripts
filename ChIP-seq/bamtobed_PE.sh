for i in `ls --color=auto ../*.bam`; do 
	tag=`echo $i | awk -F"_" '{print$2}'`;     	
	echo $tag;
	
	samtools sort -n $i $tag"_sort" ;
	bedtools bamtobed -i $tag"_sort.bam" -bedpe > $tag"_rmdup.tagAlign" ; 


done
