#input arquivo .bed ou .tagAligh referente a cada amostra ou input de cada amostra

num_r=`cat $1 | wc -l`;
name_r=`echo $1| awk -F"." '{print$1}'`;
echo "Número de reads de "$name_r": "$num_r;

perl CalBedNrf.pl $1;

if[$2]; then 
 IntersectBed -a $1 -b $2 -c -f 0.20  > Output.intersectBed; 
 perl getCnt.pl;
 rm -f Output.intersectBed;
fi


