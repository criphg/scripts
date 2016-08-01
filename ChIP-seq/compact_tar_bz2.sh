#compactação máxima usando tar e bz2

for i in `ls $1`; do 

	tar cvf temp.tar $i;
	bzip2 -9 temp.tar;

done;
