#!/bin/bash
#Programa feito para randomizar e dividir amostras paired-end de chip-seq Illumina em duas replicatas técnicas.
#Necessário definir o tamanho das replicatas na variável num_out.
#ter um arquivo list com o nome dos diretórios onde estão os fastq zipados
#ter o nome do arquivo com o padrao de formato illumina paired-end _R1_ e _R2_

#for i in `cat list`; do
i=$1

num_out="1070000";

cd $i;
	mkdir Experimental_Replicates;
	cd Experimental_Replicates;

		for a in `ls --color=auto ../*_R1_*.fastq.gz`; do

			#definindo nomes de entrada e saída
			R2_x=`echo $a | sed 's/_R1_/_R2_/g'`;
			inname_R1=`echo $a | sed 's/.fastq.gz/.fastq/g'`;
			inname_R2=`echo $R2_x | sed 's/.fastq.gz/.fastq/g'`;
			outname_R1=`echo $a | awk -F"/" '{print$NF}' | sed 's/.fastq.gz/_REP1fastq/g'`; 		
			outname_R2=`echo $R2_x | awk -F"/" '{print$NF}' | sed 's/.fastq.gz/_REP1fastq/g'`;
			outname_R1_Rep2=`echo $a | awk -F"/" '{print$NF}' | sed 's/.fastq.gz/_REP2fastq/g'`;
			outname_R2_Rep2=`echo $R2_x | awk -F"/" '{print$NF}' | sed 's/.fastq.gz/_REP2fastq/g'`;
		echo "\$a= "$a;
		echo "\$R2_x= "$R2_x;

		echo "inname_R1= "$inname_R1;
		echo "inname_R2= "$inname_R2;
		echo "outname_R1= "$outname_R1;
		echo "outname_R2= "$outname_R2;
		echo "outname_R1_Rep2 = "$outname_R1_Rep2;
		echo "outname_R2_Rep2 = "$outname_R2_Rep2;

		#descompactando arquivos sem remover os originais.
		echo "Descompactando $a";
		gunzip -k $a;

		echo "Descompactando $R2_x";
		gunzip -k $R2_x;

		echo "Shuf Rep1";
		#shuf 1 - randomizando e separando replicata 1
		cat $inname_R1 | grep "^@NF" | awk '{print$1}' | shuf -n $num_out -o list_temp1;
		

		LC_ALL=C fgrep -w -A 3 -f list_temp1 $inname_R2 > $outname_R2 &
		LC_ALL=C fgrep -w -A 3 -f list_temp1 $inname_R1 > $outname_R1; wait
		sed -i "/^--$/d" $outname_R1 &
                sed -i "/^--$/d" $outname_R2; wait
		#rm -f list_temp1;

		echo "Shuf Rep2";
		#shuf 2	- randomizando e separando replicata 2
		cat $inname_R1 | grep "^@NF" | awk '{print$1}' | shuf -n $num_out -o list_temp2;
		LC_ALL=C fgrep -w -A 3 -f list_temp2 $inname_R2 > $outname_R2_Rep2 &
                LC_ALL=C fgrep -w -A 3 -f list_temp2 $inname_R1 > $outname_R1_Rep2; wait
		sed -i "/^--$/d" $outname_R1_Rep2 &
		sed -i "/^--$/d" $outname_R2_Rep2;
		#rm -f list_temp2;

		#zipando arquivos de saída.
		gzip $outname_R1 &
		gzip $outname_R2 &
		gzip $outname_R1_Rep2 &
		gzip $outname_R2_Rep2; wait

		#removendo arquivos fastq descompactados.
#		rm -f ../*.fastq ;
 done;

	cd .. ;
cd .. ;
#done 
