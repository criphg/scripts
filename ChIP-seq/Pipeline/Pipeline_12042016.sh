#Pipeline ChIP-seq Rubens 14/07/2015 - revisado 10/02/2016
#Testado - até dia 29/01/2016 - tudo ok até MERGE.

###########IMPORTANTE###########
#NOME DOS ARQUIVOS INPUTS DEVEM SER O MESMO NOME DOS CHIP ATÉ O "_" COM INSERÇÃO DE "-I"  >> ISSO É IMPORTANTE PARA SE DEFINIR NA FUNÇÃO ENRICHMENT O NOME DO CHIP E CONTROLE.
#Ex. se CHIP se chamar SAMPLE1_S1_L001_blabla.fq, o INPUT tem que se chamar SAMPLE1-I_S1_L001_blabla.fq {alternativa, mudar a entrada tagINPUT_r abaixo}
#

#dependências 
#  1 > fastxtoolkit >: sudo apt-get install fastx-toolkit
#  2 > BWA >: sudo apt-get install bwa
#  3 > MACS2 >: https://pypi.python.org/pypi/MACS2/ (more stable version) (version of use macs2 2.1.0.20150731)
#  4 > MAnorm - diretório MAnorm deve ser criado e adicionado os arquivos com permissão de execução editados a partir do arquivo original encontrado em  >: http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/MAnorm.htm#Downloads 
#  5 > Bedtools >: sudo apt-get install bedtools
#  6 > SAMTOOLS >: sudo apt-get install samtools
#  7 > R >: sudo apt-get update; sudo apt-get install r-base-dev
#  8 > ChIPpeakAnno >: http://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html >> dependencies: library linux libssl (Ubuntu): sudo apt-get install libssl-dev, PACOTE biocunductor AnnotationHub: source ("http://bioconductor.org/biocLite.R") ; biocLite("AnnotationHub") 
#  9 > phantompeakqualtools >: https://code.google.com/archive/p/phantompeakqualtools/
# 10 > mel-ngs >: install by git clone https://github.com/mel-astar/mel-ngs.git

#NECESSÁRIO!!!!
#1 Definir região de trimagem (Trimagem).
#2 Definir numero de threads (Mapeamento)
#3 Tamanho do fragmento para MAnorm DEFINIDO PARA CADA UMA DAS AMOSTRAS PRESENE NO LOG DO MACS.
#4 Diretório contendo scripts (MAnorm << criar permissão de execução - pode ser colocado em outro diretório desde que seja referenciado abaixo, mel-ngs e phantompeakqualtools - todos devem estar em um mesmo diretório).
#5 Diretório contendo genoma de referencia.
#6 Diretório contendo o arquivo GFF.

#=========================================================================
#Variáveis a serem definidas:

##1 - TRIMAGEM
#lista de diretórios contendo fastq.gz
list_r="list_samples";
#read length:
read_trim_lgh="100";
#diretórios cotendo fastq.gz dentro da list_samples

##2 - ALIGHMENT
#genome directory and name file
genome_r="Genome/TriTrypDB-9.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta";
#diretório de saída do BAM
outBAM_r="BWA_BAM_MEM";

#número de processos
thread_r="8";

##3 - MERGE
#diretório de saída
mergedir="Merged";

##4 - Quality filter - phred
#Qualidade phred 
qualNum_r="20";

##5 - Enrichment
#diretórios de saída para o MACS
dirmacsN_r="MACS2_Narrow";
dirmacsB_r="MACS2_Broad";
#Tag usada no nome dos arquivos para identificar os inputs
tag_INPUT_r="-I";
#Genome effective size
genomeEffSize_r="17941408";

##6 - MAnorm
#local onde está o arquivos editados do MAnorm.
dirProgMAnorm="MAnorm";
#diretório de saída para o MAnorm
dirMAnorm="Result_MAnorm";
#Nome contido no arquivo fq dado a amostra 1 (geralmente o controle, não tratado, ex. Ty)
x_tag1_r="Ty";
#Nome contido no arquivo fq dado a amostra 2 (tratado, ex. TyM)
x_tag2_r="TyM";
#cuttoffs para classificação MAnorm por M value e pvalue.
#Valor de M a ser considerado para as regiões comuns. Ex. para valor de 0.5 vai ser considerado M de -0,5 a 0,5 como pico comum.
x_CF_M_Common="0.5";
#Valor de M a ser considerado para pico único em cada uma das amostras.
x_CF_M_Uniq="1";
#Valor máximo de -log10(pvalue) a ser considerado (5 eq p-value 10^-5).
x_CF_log10pvalue="5";

##7 - IGV-files
#diretório de saída para os arquivos do IGV
dirIGV_r="IGVfiles";
#diretório onde está o arquivo GFF.
codingdir_r="GFF";
#nome do arquivo GFF.
gffname_r="TriTrypDB-9.0_TcruziCLBrenerEsmeraldo-like_gene.gff";

##8 - ChipPeakAnno - busca das regiões flanqueadoras e fasta das regiões de controle putativa
#diretório de saída para as análises
dirAnalisys_r="AnnotatePeak";
#diretório onde está os scripts R_ChIP.R e joinData2.pl 
dirScripts_r="scripts";

##9 - Métricas
#diretório de saída pra as métricas
metricsdir="metrics";
#diretório do phantompeakstools
phantomdir="scripts/phantompeakqualtools";
#diretório de saída para cálculo de FRiP e NRF
outdirNRFFrip="NRF_FRiP";
#diretório de localização do mel-ngs
melngsdir="scripts";

##10 - Análise de Motivos
#diretório de saída para análise de motivos
motifdir="Motif_Analysis";
#diretório de saída para análise do HOMER
homerdir="HOMER_RESULTS";
#tamanho do fragmento a ser análisedo para busca de domínios (opção given homer usa extamente o tamanho da região que lhe é dado como pico enriquecido). Além dessa opção há a busca pelo tamanho médio encontrado pelo MACS.
#Amostra 1
homerfragsize_S1="given";
#Amostra 2
homerfragsize_S2="given";
#diretório de saída para a análise com determinado tamanho de fragmento especificado acima dentro do diretório do HOMER especificado anteriormente ($homerdir).
homeroutdir="Result_FSgiven";




#=========================================================================
#FUNÇÕES

##1 TRIMAGEM

trimag_r () {
echo "##1## TRIMAGEM";
cat $1;

for a in `cat $1`; do
cd $a;

        for i in `ls --color=auto *fastq.gz`; do
                zcat $i | fastx_trimmer -Q33 -z -l $2 -i - -o $i"_trim"$2"bp.gz";
        done;

cd .. ;

done;
}

##2 ALIGHMENT
map_r () {
echo "##2## MAPEAMENTO";

bwa index $2;

for i in `cat $1`;do
cd $i;
        mkdir $3;
        for a in `ls --color=auto *_R1_*fastq.gz_trim$5bp.gz`;do
           r2=`echo "$a"|sed 's/_R1_/_R2_/g'`;
           sam=`echo "$a"|sed 's/_R1_/_R1-2_/g'`;

           bwa aln -q 15 -t $4 -f "./"$3"/"$a"_BWA.sai" "../"$2 $a;
           bwa aln -q 15 -t $4 -f "./"$3"/"$r2"_BWA.sai" "../"$2 $r2;
           bwa sampe -a 750 -f "./"$3"/"$sam"_BWA.sam" "../"$2 "./"$3"/"$a"_BWA.sai" "./"$3"/"$r2"_BWA.sai" $a $r2;
           samtools view -bSu "./"$3"/"$sam"_BWA.sam"| samtools sort -n -o - samtools_nsort_tmp |samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp |samtools fillmd -u - "../"$2 > "./"$3"/"$sam"_BWA.sam_fixed.bam";
           samtools flagstat "./"$3"/"$sam"_BWA.sam_fixed.bam" >> "./"$3"/"$sam"_BWA.sam_fixed.bam_flgst";
        done;
##||        rm -f *.sam *.sai ;
cd ..;
done;

}

map_BWAMEN_r () {
echo "##2## MAPEAMENTO";

bwa index $2;

for i in `cat $1`;do
cd $i;
        mkdir $3;
        for a in `ls --color=auto *_R1_*fastq.gz_trim$5bp.gz`;do

	echo "##########a = $a"
        echo "##########\$2 = $2"

           r2=`echo "$a"|sed 's/_R1_/_R2_/g'`;
           sam=`echo "$a"|sed 's/_R1_/_R1-2_/g'`;


           bwa mem -M -t $4 "../"$2 $a $r2 > "./"$3"/"$sam"_BWA.sam";
           samtools view -bSu "./"$3"/"$sam"_BWA.sam"| samtools sort -n -o - samtools_nsort_tmp |samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp |samtools fillmd -u - "../"$2 > "./"$3"/"$sam"_BWA.sam_fixed.bam";
           samtools flagstat "./"$3"/"$sam"_BWA.sam_fixed.bam" >> "./"$3"/"$sam"_BWA.sam_fixed.bam_flgst";
        done;
##||        rm -f *.sam *.sai ;
cd ..;
done;

}


##3 MERGE

merge_r () {
echo "###3### Merging"

mkdir $1;
for i in `cat $2`; do
        Name=$1"/"`ls --color=auto $i/$3/*_L001_*.bam| sed 's/L001/MERGED/g' | awk -F"/" '{print $NF}'`;
        List=`ls --color=auto $i/$3/*.bam | tr '\n' ',' | sed 's/,/ /g'`;
        OutName=`echo $Name | awk -F'/' '{print $NF}' | awk -F"_" '{print$1"_MERGED_"$4"_"$6"_sort_index"}'`;

        echo "Name = $Name";
        echo "List = $List";
        echo "OutName = $OutName";

        samtools merge $Name $List;
        samtools sort $Name $OutName;
        samtools index $OutName".bam";
done;

##|| rm -f ./Merged/*fixed.bam;

echo "==========================FIM MERGE!================================";

}


##4 QUALITY FILTER - Q30

quality_pherd_r() {

echo "###4### Filter";

cd $2;

        mkdir "Q"$1;
        for i in `ls --color=auto *index.bam`; do

                Name="Q$1/Q$1_"$i"_tempX"
                Name2=`a="Q$1/Q$1_"$i ; echo ${a%.*m}`;
        #       echo "Name =" $Name;
        #       echo "Name 2 = "$Name2;
                samtools view -q$1 -b $i > $Name;
                samtools sort $Name $Name2;
                samtools index $Name2".bam";
                samtools flagstat $Name2".bam" > $Name2"_flgst";
        done;

        rm -f Q$1/*_tempX;

cd .. ;

}


##4.b QUALITY FILTER PROPER PAIR AND BOTH PAIRS UNIQUE MAPPED.
################FALTA EDITAR!!!!! COLOCAR JUNTO COM O DE CIMA E EDITAR ENRICHMENT
#Acreditei no flagstat... possível problema com a seleção dos pares... verificando.
# Filtro dos arquivos de  qualidade Q30, por propriamente mapeados -f 3, onde ambos os membros do par são unicamente mapeados.

filter_mateBWA_XTAR_r() {

echo "###4.b### QUALITY PAIR";

cd $2"/Q"$1 ;

mkdir "Q"$1"_pMapped";

        for i in `ls --color=auto "Q$1_*index.bam"`; do

        cd "Q"$1"_pMapped";
                outname=`echo $i | awk -F"." '{print$1}'`;
                echo "Outname =" $outname;
                samtools view -h -f 3 "../"$i > "temp_file" ;
                cat temp_file | grep -v "^@" | awk '{if ($12=="XT:A:R") print$1}' > list_temp;
                cat temp_file | fgrep -v -f list_temp > temp_file2;
                samtools view -bS temp_file2 > temp_file3
                samtools sort temp_file3 $outname"_pMapped_index";
                samtools index $outname"_pMapped_index.bam";
                samtools flagstat $outname"_pMapped_index.bam" > $outname"_pMapped_flgst";

                rm -f temp_file* ;
                rm -f list_temp;
        cd ..;
        done;

cd ../../;

}

##5 Enrichment

##NARROW
x_macs2Narrow_r() {

echo "##5## Enrichment_NARROW"


mkdir $1;

for i in `ls --color=auto ${3}/Q${2}/Q${2}_pMapped/Q${2}_*index.bam | grep '\\'${4}`; do
	echo $i;
        control=`echo "../$i"`;
        chip=`echo "../"$i | sed 's/\-I//g'`;
        name=`echo $chip | awk -F"/" '{print$NF}' | awk -F"." '{print$1}'`;

#name=`echo $name | awk -F"_" '{print$2}'`;

        echo "Control==="$i
        echo "ChIP ==="$chip;
        echo "Name === "$name;

        #Acrescentar -B para obtenção do bedgraph, não é possível geração do bedgraph e dos picos na mesma linha análise.
        #Para outros instalações pode ser necessário mudar o nome macs para macs14.
        #cd MACS14_Q30 ;                
        #       macs -S -w -t $chip -c $control -g 32529070 -n $name -f BAM --call-subpeaks --verbose 3 &> $name"_macs14_wig.log";
        #       macs -S -B -t $chip -c $control -g 32529070 -n $name -f BAM --call-subpeaks --verbose 3 &> $name"_macs14_bdg.log";
        #cd .. ;


        ###Narrow
        #cd MACS2_N_dup/ ;
        #      macs2 callpeak -g 25066458 -n $name -B --keep-dup 'all' --call-summits -f BAMPE -t $chip -c $control &> $name"_n.log";
        #cd .. ;

        ###Narrow _ NO summits _ subdivisão de regiões enriquecidas gera problemas na análise posterior usando MAnorm, assim optou-se por não usar --call-summits.

        cd $1/ ;
                macs2 callpeak -g $5 -n $name -B --keep-dup 'all' -f BAMPE -t $chip -c $control &> $name"_n.log";
        cd .. ;

done;

}

##BROAD
x_macs2Broad_r() {

echo "##5## Enrichment_BROAD"


mkdir $1;

for i in `ls --color=auto ${3}/Q${2}/Q${2}_pMapped/Q${2}_*index.bam | grep '\\'${4}`; do
        control=`echo "../$i"`;
        chip=`echo "../"$i | sed 's/\-I//g'`;
        name=`echo $chip | awk -F"/" '{print$NF}' | awk -F"." '{print$1}'`;

#name=`echo $name | awk -F"_" '{print$2}'`;

        echo "Control==="$control;
        echo "ChIP ==="$chip;
        echo "Name === "$name;




###Broad

        cd $1/ ;
                macs2 callpeak -g $5 --keep-dup 'all' -n $name -t $chip -c $control -f BAMPE --broad -B &> $name"_b.log";
        cd .. ;

done;

}


##6 - MAnorm

# Tranformar BAMtoBED para busca de reads.

BAMtoBED_MA_reads_r() {
echo "##Getting Reads BAMtoBED##"
echo"";

cd "${2}/Q${1}/Q${1}_pMapped/";

        mkdir bed;

        for i in `ls --color=auto Q${1}_*pMapped_index.bam | grep -v '\\'${3}`; do
                arq_name=`echo $i | awk -F"." '{print"bed/"$1".bed"}'`;
                echo "Name ="$arq_name;

                bedtools bamtobed -i $i > $arq_name

        done;


        #Formatando arquivos .bed para o formato necessário
	echo "Formatando arquivos para formato BED para MAnorm...";
        cd bed;



                for i in `ls --color=auto Q${1}_*.bed`; do


###
### <<<<<<<<<<<<<<<<EDITAR TAG ABAIXO!!!!!!
###
                        tag=`echo $i | awk -F "_" '{print$2}'`
			echo "Tag file="$tag
                        awk '{print$1"\t"$2"\t"$3"\t"$6}' $i > "reads_"$tag".bed";
         #              awk '{gsub("TcC","c",$1); gsub("-S","",$1); print$1"\t"$2"\t"$3"\t"$6}' $i > "reads_"$tag".bed";

                done;
        cd ..;
cd ../../../;
}

BAMtoBED_MA_reads_INPUTS_r() {
echo "##Getting Reads BAMtoBED_INPUT##"
echo"";

cd "${2}/Q${1}/Q${1}_pMapped/";

        mkdir bed;
	mkdir bed/INPUT;
	
        for i in `ls --color=auto Q${1}_*pMapped_index.bam | grep '\\'${3}`; do
                arq_name=`echo $i | awk -F"." '{print"bed/INPUT/"$1".bed"}'`;
                echo "Name ="$arq_name;

                bedtools bamtobed -i $i > $arq_name

        done;


        #Formatando arquivos .bed para o formato necessário
        echo "Formatando arquivos para formato BED para MAnorm...";
        cd bed/INPUT;



                for i in `ls --color=auto Q${1}_*.bed`; do


###
### <<<<<<<<<<<<<<<<EDITAR TAG ABAIXO!!!!!!
###
                        tag=`echo $i | awk -F "_" '{print$2}'`
                        echo "Tag file="$tag
                        awk '{print$1"\t"$2"\t"$3"\t"$6}' $i > "reads_"$tag".bed";
         #              awk '{gsub("TcC","c",$1); gsub("-S","",$1); print$1"\t"$2"\t"$3"\t"$6}' $i > "reads_"$tag".bed";

                done;
        cd ../../;
cd ../../../;
}








# Formatar Peaks no MACS

MA_peaks_r(){


echo "##Getting Peaks to MAnorm##";
echo "";

for i in `ls --color=auto {${1},${2}}/Q${3}_*peaks* | grep "broadPeak\|narrowPeak"`; do

	echo "MACS file="$i ;
        dir=`echo $i | awk -F "/" '{print$1}'`;

        cd $dir;
                tag=`echo $i | awk -F "/" '{print$NF}' | awk -F[_.]+ '{print$2"_"$NF}' | sed -e 's/broadPeak/B/g;s/narrowPeak/N/g'`;
                file=`echo $i | awk -F "/" '{print$NF}'`;
                echo "Peak_tag="$tag;
                echo "Peak_file="$file
        #       awk '{gsub("TcC","c",$1); gsub("-S","",$1); print$1"\t"$2"\t"$3}' $file > "peaks_"$tag;
                awk '{print$1"\t"$2"\t"$3}' $file > "peaks_"$tag".bed";

        cd .. ;

done;

}



##Execução do MAnorm

MAnorm_r(){


echo "###6### MAnorm";
echo "";

######DEVE HAVER DIRETÓRIO MAnorm CRIADO COM OS ARQUIVOS MAnorm.sh (editado) e MAnorm.r com permissão de execução nele.

##Busca de tamanho de fragmentos:
#Tamanho do fragmento deve ser recuperado antes de fazer MAnorm (está descrito no log do MACS na frente do fragment size = , linhas abaixo recuperam esse valor para cada amostra.


fragBroad_Sample1=`cat ${4}/Q${1}_${7}_*_index_b.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
fragBroad_Sample2=`cat ${4}/Q${1}_${8}_*_index_b.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
fragNarrow_Sample1=`cat ${3}/Q${1}_${7}_*_index_n.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
fragNarrow_Sample2=`cat ${3}/Q${1}_${8}_*_index_n.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;

echo "##Fragment size##";
echo "${4}/Q${1}_${7}_*_index_b.log";
echo "$fragBroad_Sample1";
echo "${4}/Q${1}_${8}_*_index_b.log";
echo "$fragBroad_Sample2";
echo "${3}/Q${1}_${7}_*_index_n.log";
echo "$fragNarrow_Sample1";
echo "${3}/Q${1}_${8}_*_index_n.log";
echo "$fragNarrow_Sample2";
echo""

mkdir ${5};

cd ${5};

        cp ../${6}/{MAnorm.sh,MAnorm.r,classfy_by_M_value.sh} .; #MAnorm é o diretório de origem do programa.

        ln -s ../${2}/Q${1}/Q${1}_pMapped/bed/reads_*.bed . ;
        ln -s ../{${3},${4}}/peaks_*.bed . ;

        mkdir Broad;
        mkdir Narrow;

        cd Broad;
                cp ../MAnorm.r . ;
                cp ../MAnorm.sh . ;
                cp ../classfy_by_M_value.sh . ;

                ./MAnorm.sh ../peaks_${7}_B.bed ../peaks_${8}_B.bed ../reads_${7}.bed ../reads_${8}.bed $fragBroad_Sample1 $fragBroad_Sample2 ;
                ./classfy_by_M_value.sh ${9} ${10} ${11};
        cd ..;



        cd Narrow;
                cp ../MAnorm.r . ;
                cp ../MAnorm.sh . ;
                cp ../classfy_by_M_value.sh . ;

                ./MAnorm.sh ../peaks_${7}_N.bed ../peaks_${8}_N.bed ../reads_${7}.bed ../reads_${8}.bed $fragNarrow_Sample1 $fragNarrow_Sample2 ;
                ./classfy_by_M_value.sh ${9} ${10} ${11};
        cd ..;
cd ..;

}


##ARQUIVO GFF PRECISA SER PREPARADO. ELIMINAÇÃO DO CABEÇALHO E SELEÇÃO SOMENTE DOS GENES.

format_gff_r(){
cd ${1};

	grep -v "^#" ${2} | awk '{if($3 == "gene") print$0}' > Genes.gff;

cd .. ;
}




igv_files_r(){



############## IGV-files

echo "###7### IGV-files";

echo "##Copiando e renomeando arquivos para o diretório $1 ##";
mkdir ${1};


cd ${1};


#########EDITAR ARQUIVOS PARA ADMIIR QUALQUER TAG

        mkdir Broad;
        mkdir Narrow;

#echo "###########111";
#broad
        cp ../${2}/Broad/{sample*.bed,unbiased*.bed,*.wig} ./Broad;
        mv Broad/sample1_uniq_peaks.bed Broad/Uniq${3}_MAnorm_B.bed;
        mv Broad/sample2_uniq_peaks.bed Broad/Uniq${4}_MAnorm_B.bed;
        mv Broad/sample1_peaks.wig Broad/${3}_MAnorm_Mval_B.wig;
        mv Broad/sample2_peaks.wig Broad/${4}_MAnorm_Mval_B.wig;
        mv Broad/unbiased_peaks.bed Broad/CommonPeak_B.bed;

#echo "#############222";
        cp ../${6}/{*.bdg,*.broadPeak} ./Broad;
        mv ./Broad/Q${7}_${3}_MERGED*_pileup.bdg ./Broad/ChIP_${3}_MACS_B.bdg;
        mv ./Broad/Q${7}_${4}_MERGED*_pileup.bdg ./Broad/ChIP_${4}_MACS_B.bdg;
        mv ./Broad/Q${7}_${3}_MERGED*_lambda.bdg ./Broad/INPUT_${3}_MACS_B.bdg;
        mv ./Broad/Q${7}_${4}_MERGED*_lambda.bdg ./Broad/INPUT_${4}_MACS_B.bdg;

        mv ./Broad/Q${7}_${3}_MERGED*.broadPeak ./Broad/${3}_MACS_B.broadPeak;
        mv ./Broad/Q${7}_${4}_MERGED*.broadPeak ./Broad/${4}_MACS_B.broadPeak;


        cp ../${9}/Genes.gff ./Broad;

#echo "#############333"

#narrow
        cp ../${2}/Narrow/{sample*.bed,unbiased*.bed,*.wig} ./Narrow;
        mv ./Narrow/sample1_uniq_peaks.bed ./Narrow/Uniq${3}_MAnorm_N.bed;
        mv ./Narrow/sample2_uniq_peaks.bed ./Narrow/Uniq${4}_MAnorm_N.bed;
	mv ./Narrow/sample1_peaks.wig ./Narrow/${3}_MAnorm_Mval_N.wig;
        mv ./Narrow/sample2_peaks.wig ./Narrow/${4}_MAnorm_Mval_N.wig;
        mv ./Narrow/unbiased_peaks.bed ./Narrow/CommonPeak_N.bed;

        cp ../${5}/{*.bdg,*.narrowPeak} ./Narrow;
        mv ./Narrow/Q${7}_${3}_MERGED*_pileup.bdg ./Narrow/ChIP_${3}_MACS_N.bdg;
        mv ./Narrow/Q${7}_${4}_MERGED*_pileup.bdg ./Narrow/ChIP_${4}_MACS_N.bdg;
        mv ./Narrow/Q${7}_${3}_MERGED*_lambda.bdg ./Narrow/INPUT_${3}_MACS_N.bdg;
        mv ./Narrow/Q${7}_${4}_MERGED*_lambda.bdg ./Narrow/INPUT_${4}_MACS_N.bdg;

        mv ./Narrow/Q${7}_${3}_MERGED*.narrowPeak ./Narrow/${3}_MACS_N.narrowPeak;
        mv ./Narrow/Q${7}_${4}_MERGED*.narrowPeak ./Narrow/${4}_MACS_N.narrowPeak;
        cp ../${9}/Genes.gff ./Narrow;

cd ..;
}


#Função acessória de x_peakanalisys_r para formatar os arquivos de saída.

clean_r () {
	awk -F"\t" -v col_=18 '{ 

        gsub("[+]", " ",$col_);
        gsub("%2F", "/",$col_);
        gsub("%5B", "[",$col_);
        gsub("%5D", "]",$col_);
        gsub("%27", "`",$col_);
        gsub("%2C", ",",$col_);
        gsub("%28", "(",$col_);
        gsub("%29", ")",$col_);
        gsub("description=", "",$col_);
        gsub("##Name=", "",$16);
        OFS="\t";
        print$0;

        }' $1 > foo_xx

	tag=`echo $1 | awk -F"_" '{print$1}'`;
  	cat foo_xx | head -n 1 | awk '{$2="Chr" FS $2}1' | awk '{gsub("\"", ""); print$0}' | awk '{$16="Analisys_Type"; $17="Gene_ID"; $18="Closest_gene_name"; $19="";print$0}' OFS="\t" > "foo_ED";

	tail -n +2 foo_xx | awk -F"\t" 'BEGIN{OFS="\t"}{gsub("\"","");gsub(" ","_",$6); print$0}' >> "foo_ED";

	cat foo_ED > $1

  rm -rf foo_xx foo_ED;
}





#função acessória de x_peakanalisys_r para formatar os arquivos de saída. 
agrup_r(){
	#$1 = ChIPpeakAnnoFile *.xls
	head -1 $1 > head1;
	#$2 = MAnormfile Produzido pelo classify editado... S1_UNIQ_PEAK.out ou S2_UNIQ_PEAK.out
	head -1 $2 > head2;
	paste -d"\t" head2 head1 > head;

	sed '1d' $1 | sort -k8 -V > temp1
	sed '1d' $2 | sort -k1 -V > temp2
	paste -d"\t" temp2 temp1 > temp3;

	#$3 = nome do arquivo de saída | Usei cut para remover a coluna space.
	cat head temp3 | cut -f-10,12- > $3;
	
	rm -f temp1 temp2 temp3;
}


x_peakanalisys_r(){





################5 ChIPeakAnno

echo "###8### Peak analysis";

echo "ChIPeakAnno...";
echo "Broad...";

mkdir ${1};

cd ${1};

        mkdir Broad;
        mkdir Narrow;

        cd Broad;
                cp ../../${2}/R_ChIP.R .;
                cp ../../${2}/joinData2.pl .;

                Rscript R_ChIP.R ../../${3}/Genes.gff "../../"${5}"/Broad/sample1_uniq_peaks.bed" ${6}"_Broad";
                perl joinData2.pl ../../${3}/Genes.gff ${6}"_Broad.aPeaks" > ${6}"_Broad_aP_join.xls";

                Rscript R_ChIP.R ../../${3}/Genes.gff "../../"${5}"/Broad/sample2_uniq_peaks.bed" ${7}"_Broad";
                perl joinData2.pl ../../${3}/Genes.gff ${7}"_Broad.aPeaks" > ${7}"_Broad_aP_join.xls";

                #Limpando xls e juntando com arquivo contendo informa
                clean_r ${6}"_Broad_aP_join.xls";
                clean_r ${7}"_Broad_aP_join.xls";


		#Agrupando arquivos CHIPPEAKANNO com Manorm


		agrup_r ${6}"_Broad_aP_join.xls" "../../"${5}"/Broad/S1_UNIQ_PEAK.out" ${6}"_Broad_tabFull.csv";
		agrup_r ${7}"_Broad_aP_join.xls" "../../"${5}"/Broad/S2_UNIQ_PEAK.out" ${7}"_Broad_tabFull.csv";







echo "Fasta file...";
                mkdir Fasta;
                cd Fasta;

                        cat "../"${6}"_Broad.aPeaks" | tr -d "\"" | awk '{print $2":"$3"-"$4}' | tail -n +2 > TEMP_FILE.faidx;
                        for i in `cat TEMP_FILE.faidx`; do samtools faidx ../../../${8} $i; done > TEMP_FILE_2.faidx;


#Remove quebra de linhas das sequencias fasta (facilita a feitura do primer)
                        awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' TEMP_FILE_2.faidx > ${6}"_Broad_aP.fasta";

                        rm -f {TEMP_FILE.faidx,TEMP_FILE_2.faidx};

                        cat "../"${7}"_Broad.aPeaks" | tr -d "\"" | awk '{print $2":"$3"-"$4}' | tail -n +2 > TEMP_FILE.faidx
                        for i in `cat TEMP_FILE.faidx`; do samtools faidx ../../../${8} $i; done > TEMP_FILE_2.faidx;
                        awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' TEMP_FILE_2.faidx > ${7}"_Broad_aP.fasta";

			rm -f {TEMP_FILE.faidx,TEMP_FILE_2.faidx};

                cd .. ;

        cd .. ;


#########################################################################
echo "Narrow...";
        cd Narrow;

                cp ../../${2}/R_ChIP.R .;
                cp ../../${2}/joinData2.pl .;

                Rscript R_ChIP.R ../../${3}/Genes.gff "../../"${5}"/Narrow/sample1_uniq_peaks.bed" ${6}"_Narrow";
                perl joinData2.pl ../../${3}/Genes.gff ${6}"_Narrow.aPeaks" > ${6}"_Narrow_aP_join.xls";

                Rscript R_ChIP.R ../../${3}/Genes.gff "../../"${5}"/Narrow/sample2_uniq_peaks.bed" ${7}"_Narrow";
                perl joinData2.pl ../../${3}/Genes.gff ${7}"_Narrow.aPeaks" > ${7}"_Narrow_aP_join.xls";


		#Limpando xls
		clean_r ${6}"_Narrow_aP_join.xls";
		clean_r ${7}"_Narrow_aP_join.xls";


                #Agrupando arquivos CHIPPEAKANNO com Manorm


                agrup_r ${6}"_Narrow_aP_join.xls" "../../"${5}"/Narrow/S1_UNIQ_PEAK.out" ${6}"_Narrow_tabFull.csv";
                agrup_r ${7}"_Narrow_aP_join.xls" "../../"${5}"/Narrow/S2_UNIQ_PEAK.out" ${7}"_Narrow_tabFull.csv";






echo "Fasta file...";

                mkdir Fasta;
                cd Fasta;

                        cat "../"${6}"_Narrow.aPeaks" | tr -d "\"" | awk '{print $2":"$3"-"$4}' | tail -n +2 > TEMP_FILE.faidx;
                        for i in `cat TEMP_FILE.faidx`; do samtools faidx ../../../${8} $i; done > TEMP_FILE_2.faidx;
                        awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' TEMP_FILE_2.faidx > ${6}"_Narrow_aP.fasta";

                        rm -f {TEMP_FILE.faidx,TEMP_FILE_2.faidx};

                        cat "../"${7}"_Narrow.aPeaks" | tr -d "\"" | awk '{print $2":"$3"-"$4}' | tail -n +2 > TEMP_FILE.faidx;
                        for i in `cat TEMP_FILE.faidx`; do samtools faidx ../../../${8} $i; done > TEMP_FILE_2.faidx;
                        awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' TEMP_FILE_2.faidx > ${7}"_Narrow_aP.fasta";

                        rm -f {TEMP_FILE.faidx,TEMP_FILE_2.faidx};
                cd .. ;
        cd .. ;
cd ..;
}

##8b - incorporar funcoes para formataćão de tabela - 1- comparar arquivos joinData.pl 1 e 2, 2 - acrescentar funções localizadas no arquivo :/media/sf_D_DRIVE/Dados/Tcruzi/Novos_primers_20012016$ vim Paste_MAnornResultsGFF_22012016.sh







##9 - Métricas
#formatar arquivo SPP localizado na pasta /media/sf_D_DRIVE/Dados/Tcruzi/SPP_Analise , arquivo SPP_quality_analisys.R
#verificar uso do runSPP.R para calculo dessa métrica /Tools/phantompeakqualtools
#Fazer o link simbólico do Ty.bed para Ty.tagAlign
#verificar outras métricas em ~/Tools/mel-ngs/mel-chipseq/chipseq-metrics 


#Não é necessário separar para broad e narrow, SPP leve em consideração somente reads mapeadas para buscar Cross-correlation.
x_SPP_r () {

#cria diretório de métricas, se já existir não mostra erro.
mkdir -p ${4};
cd ${4};


	echo "##9## Metrics - SPP - Cross-correlation";

	#Busca dos arquivos bed, link simbólico mudando o nome do arquivo para .tagAlign >> necessário para o SPP
	for i in `ls ../${2}/Q${1}/Q${1}_pMapped/{bed,bed/INPUT}/read*.bed`; do
        	namef=`echo $i | awk -F"/" '{print$NF}' | awk -F"." '{print$1}' | awk -F"reads_" '{print$2}'`;
        	echo $namef;
        	ln -s $i ./${namef}".tagAlign" ;
	done;


	mkdir SPP_CrossCorrelation;
	cd SPP_CrossCorrelation;

        	#Busca SPP.R

        	ln -s ../../${3}/run_spp.R .;


		#Arquivo de Explicação da saída

		echo "
#===========================
#TYPICAL USAGE
#===========================
#(1) Determine strand cross-correlation peak / predominant fragment length OR print out quality measures
#       
#       Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>
#
#-savp will create a pdf showing the cross-correlation plot
#-out=<outFile> will create and/or append to a file named <outFile> several important characteristics of the dataset.
#The file contains 11 tab delimited columns
#
#COL1: Filename: tagAlign/BAM filename
#COL2: numReads: effective sequencing depth i.e. total number of mapped reads in input file
#COL3: estFragLen: comma separated strand cross-correlation peak(s) in decreasing order of correlation.
#         The top 3 local maxima locations that are within 90% of the maximum cross-correlation value are output.
#         In almost all cases, the top (first) value in the list represents the predominant fragment length.
#         If you want to keep only the top value simply run
#         sed -r 's/,[^\t]+//g' <outFile> > <newOutFile>
#COL4: corr_estFragLen: comma separated strand cross-correlation value(s) in decreasing order (col2 follows the same order)
#COL5: phantomPeak: Read length/phantom peak strand shift
#COL6: corr_phantomPeak: Correlation value at phantom peak
#COL7: argmin_corr: strand shift at which cross-correlation is lowest
#COL8: min_corr: minimum value of cross-correlation
#COL9: Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8
#COL10: Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8)
#COL11: QualityTag: Quality tag based on thresholded RSC (codes: -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh)
#
#You can run the program on multiple datasets in parallel and append all the quality information to the same <outFile> for a summary analysis.
#
#NSC values range from a minimum of 1 to larger positive numbers. 1.1 is the critical threshold. 
#Datasets with NSC values much less than 1.1 (< 1.05) tend to have low signal to noise or few peaks (this could be biological eg.a factor that truly binds only a few sites in a particular tissue type OR it
# could be due to poor quality)
#
#RSC values range from 0 to larger positive values. 1 is the critical threshold.
#RSC values significantly lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mismappings, shallow 
#sequencing depth (significantly below saturation) or a combination of these. Like the NSC, datasets with few binding sites (< 200) which is biologically justifiable also show low RSC scores.
#
#Qtag is a thresholded version of RSC.
############Ref: README.txt, phantompeakqualtools, https://code.google.com/archive/p/phantompeakqualtools/ 
                        
" > "README_BasicExplanationSPPOutput";


        	#Usando SPP (resultados são os mesmos usando ou não o INPUT para cálculo de cross-correlação de ChIP, por isso não usei o INPUT e calculei sua cross-correlação separadamente).
        	for a in `ls --color=auto ../*.tagAlign`; do

                	namefx=`echo $a | awk -F"/" {'print$NF'} | awk -F"." '{print$1}'`;
                	echo "$a";
                	echo "namefx = "$namefx;

			Rscript run_spp.R -c=$a -savp -out=temp_file &> log_SPP_$namefx;

                        #cabeçalho do arquivo de saída.
                        echo -e "File_Name\tReads_Num\testFragLen\tcorr_estFragLen\tphantomPeak\tcorr_phantomPeak\targmin_corr\tmin_corr\tNSC\tRSC\tQualityTag" >> $namefx"_SPP.out";
                        cat temp_file >> $namefx"_SPP.out";

               rm -f temp_file;

                done;

	mv ../*.pdf .;
        cd .. ;
cd ..;
}


x_NRF_r() {

echo "###9### Metrics - NRF";
mkdir -p ${1};
cd ${1};

	mkdir NRF;
	cd NRF;

		echo -e "File_NAME\tNum_of_pos_uniqMapped\tuniqMapped\tNRF" > results_NRF.out;
		for b in `ls --color=auto ../*.tagAlign`; do
		
			num_r=`cat $b | wc -l`;
			name_r=`echo $b| awk -F"/" '{print$NF}' | awk -F"." '{print$1}'`;
			echo "Número de reads de "$name_r": "$num_r;
		##NRF
			perl ../../${2}/mel-ngs/mel-chipseq/chipseq-metrics/CalBedNrf.pl $b >> results_NRF.out;
		done;
	cd ..;
cd .. ;

}

x_FRiP_N_r(){

echo "###9### Metrics - FRiP_N";


mkdir -p ${1};
cd ${1};

        mkdir -p FRiP;
        cd FRiP;

		echo -e "File_Name\tUniq_reads\tPeakInside_reads\tFRiP" > FRiP_result_Narrow.tab; 
                for b in `ls --color=auto ../*.tagAlign`; do

                        num_r=`cat $b | wc -l`;
                        name_r=`echo $b| awk -F"/" '{print$NF}' | awk -F"." '{print$1}'`;
                        echo "Número de reads de "$name_r": "$num_r;



                ##FRiP

                        #narrowa
			#a variável ajust é para usar os picos de Ty e TyM para cálculo do FRiP de Ty-I e TyM-I também (interessante verificar a fração de reads do INPUT que cai dentro dos picos identificados, espera-se que seja menor do que a fração das reads_ChIP).
			ajust=`echo $name_r | sed "s/${4}//g"`;
                        peakname_N=`ls ../../${3}/*${ajust}_*.narrowPeak`;
                        echo "peakname_N = "$peakname_N;
			
			echo "intersecbed"
			echo "b= "$b;
			echo "peakname= "$peakname_N;
                        intersectBed -a $b -b $peakname_N -c -f 0.20  > ${name_r}".intersectBed";
                        var=`perl ../../${2}/mel-ngs/mel-chipseq/chipseq-metrics/getCnt.pl ${name_r}.intersectBed`;
			var2=`echo $var | awk '{gsub(" ","",$2); print$2}'`;
			echo "Var= "$var2;
			#Cálculo de FRiP - também pode ser feito com python:  a=`python -c "print 20+5/2.0"` 
			#FrIP_r=`bc <<< 'scale=2; ${num_r} / ${var2}'`;
			FrIP_r=`python -c "print round(float(($var2*1.00)/$num_r),4)"`;
			#
			echo -e "Frip= "$FrIP_r;
			echo -e $b"\t"$num_r"\t"$var2"\t"$FrIP_r >> FRiP_result_Narrow.tab;
#                        rm -f Output.intersectBed;



                done;
        cd ..;

cd .. ;



}


x_FRiP_B_r(){

echo "###9### Metrics - FRiP_B";


mkdir -p ${1};
cd ${1};

        mkdir -p FRiP;
        cd FRiP;

                echo -e "File_Name\tUniq_reads\tPeakInside_reads\tFRiP" > FRiP_result_Broad.tab;
                for b in `ls --color=auto ../*.tagAlign`; do

                        num_r=`cat $b | wc -l`;
                        name_r=`echo $b| awk -F"/" '{print$NF}' | awk -F"." '{print$1}'`;
                        echo "Número de reads de "$name_r": "$num_r;



                ##FRiP
                        #usando picos totais do MACS relativo a cada amostra
                        #broad
 			#a variável ajust é para usar os picos de Ty e TyM para cálculo do FRiP de Ty-I e TyM-I também (interessante verificar a fração de reads do INPUT que cai dentro dos picos identificados, espera-se que seja menor do que a fração das reads_ChIP).
                        ajust=`echo $name_r | sed "s/${4}//g"`;
                        peakname_B=`ls ../../${3}/*${ajust}_*.broadPeak`;
                        echo "peakname_B = "$peakname_B;
		
                        intersectBed -a $b -b $peakname_B -c -f 0.20  > ${name_r}".intersectBed";
                        var=`perl ../../${2}/mel-ngs/mel-chipseq/chipseq-metrics/getCnt.pl ${name_r}.intersectBed`;
                        var2=`echo $var | awk '{gsub(" ","",$2); print$2}'`;
                        echo "Var= "$var2;
                        #Cálculo de FRiP - também pode ser feito com python:  a=`python -c "print 20+5/2.0"` 
                        #FrIP_r=`bc <<< 'scale=2; ${num_r} / ${var2}'`;
                        FrIP_r=`python -c "print round(float(($var2*1.00)/$num_r),4)"`;
                        #
                        echo -e "Frip= "$FrIP_r;
                        echo -e $b"\t"$num_r"\t"$var2"\t"$FrIP_r >> FRiP_result_Broad.tab;
#                        rm -f Output.intersectBed;




                done;
        cd ..;

cd .. ;



}




##10 - Busca de motivos 



##10a - função HOMER



x_homer_r () {

echo "Motif searching - HOMER";
mkdir -p $1;
cd $1;	
	#transformar o arquivo do MAnorm S1_UNIQ ou S2_UNIQ em bed
	cat ../../$2/$3/S1_UNIQ_PEAK.out | awk -F"\t" 'BEGIN{OFS="\t"} {print$2,$3,$4,$1}' > S1_UNIQ_PEAK_${3}.bed;
        cat ../../$2/$3/S2_UNIQ_PEAK.out | awk -F"\t" 'BEGIN{OFS="\t"} {print$2,$3,$4,$1}' > S2_UNIQ_PEAK_${3}.bed;

	#limha de comando para homer
	#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
	findMotifsGenome.pl S1_UNIQ_PEAK_${3}.bed "../../"$7 $3'/'$6'/'$8 -size $4
	findMotifsGenome.pl S2_UNIQ_PEAK_${3}.bed "../../"$7 $3'/'$6'/'$9 -size $5

	cat S1_UNIQ_PEAK_${3}.bed S2_UNIQ_PEAK_${3}.bed > S_ALL_UNIQ_PEAK_${3}.bed
	findMotifsGenome.pl S_ALL_UNIQ_PEAK_${3}.bed "../../"$7 $3'/'$6'/ALL' -size $5

cd ..;

}


motiffinder () {
echo "###10### - Motif Searching";

mkdir -p $motifdir;
cd $motifdir;
        #busca por domínio com parametros definidos para Broad e Narrow.
        x_homer_r $homerdir $dirMAnorm Narrow $homerfragsize_S1 $homerfragsize_S2 $homeroutdir $genome_r $x_tag1_r $x_tag2_r;
        x_homer_r $homerdir $dirMAnorm Broad $homerfragsize_S1 $homerfragsize_S2 $homeroutdir $genome_r $x_tag1_r $x_tag2_r;



	#busca por domínio com size médio definido pelo MACS para Broad e Narrow.
	#size médio MACS

	fragBroad_Sample1=`cat ../${3}/Q${1}_${4}_*_index_b.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
	fragBroad_Sample2=`cat ../${3}/Q${1}_${5}_*_index_b.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
	fragNarrow_Sample1=`cat ../${2}/Q${1}_${4}_*_index_n.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
	fragNarrow_Sample2=`cat ../${2}/Q${1}_${5}_*_index_n.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;


#$qualNum_r $dirmacsN_r $dirmacsB_r $x_tag1_r $x_tag2_r


	x_homer_r $homerdir $dirMAnorm Narrow $fragNarrow_Sample1 $fragNarrow_Sample2 Result_HOMER_FSizeMACS $genome_r $x_tag1_r $x_tag2_r;
        x_homer_r $homerdir $dirMAnorm Broad $fragBroad_Sample1 $fragBroad_Sample2 Result_HOMER_FSizeMACS $genome_r $x_tag1_r $x_tag2_r;




cd .. ;

}


#=============================================================================
#CHAMADA DE FUNÇÕES:

##1 - TRIMAGEM

#deve ser informado: 
# 1 - list_samples com os diretórios de cada amostra contendo os fastq.gz!
# 2 - tamanho da região a ser deixada.
##|| trimag_r $list_r $read_trim_lgh;


##2 - ALIGHMENT

#deve ser informado:
# 1 - list_samples com os diretórios de cada amostra contendo os fastq.gz!
# 2 - diretório do genoma e nome do arquivo.
# 3 - diretório de saída do bam

#map_r $list_r $genome_r $outBAM_r $thread_r $read_trim_lgh;
##|| map_BWAMEN_r $list_r $genome_r $outBAM_r $thread_r $read_trim_lgh;


##3 - MERGE
#deve ser informado:
# 1 - ditetório de saida para o Merge
##|| merge_r $mergedir $list_r $outBAM_r;

##4 - Quality Filer Phred
#deve ser informado:
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
##|| quality_pherd_r $qualNum_r $mergedir;

##4b - Filter mates with r1 XT:A:R, r2 XT:A:U
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
##|| filter_mateBWA_XTAR_r $qualNum_r $mergedir;

##5 - Enrichment
##NARROW
# 1 - diretório de saída para o MACS
# 2 - número a ser considerado para filtro por qualidade phred.
# 3 - diretório Merge
# 4 - Tag usada nos arquivos de INPUT para o MACS diferenciar INPUTS de CHIP.
# 5 - Genome effective size 
##|| x_macs2Narrow_r $dirmacsN_r $qualNum_r $mergedir $tag_INPUT_r $genomeEffSize_r ;

##BROAD
# 1 - diretório de saída para o MACS
# 2 - número a ser considerado para filtro por qualidade phred.
# 3 - diretório Merge
# 4 - Tag usada nos arquivos de INPUT para o MACS diferenciar INPUTS de CHIP.
# 5 - Genome effective size 
##|| x_macs2Broad_r $dirmacsB_r $qualNum_r $mergedir $tag_INPUT_r $genomeEffSize_r ;

##6 - MAnorm

##6a - Adequação do arquivo contendo reads de Ty e TyM para input do MAnorm
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
##|| BAMtoBED_MA_reads_r $qualNum_r $mergedir $tag_INPUT_r;

##6a.2 - Obtenção dos beds também para os INPUTS. Necessário para cálculos de métricas ou mesmo para vizualização no IGV.
##|| BAMtoBED_MA_reads_INPUTS_r $qualNum_r $mergedir $tag_INPUT_r;

##6b -Adequação do arquivo contendo peaks para input do MAnorm
# 1 - diretório de saída para o MACS Narrow
# 2 - diretório de saída para o MACS Broad
# 3 - número a ser considerado para filtro por qualidade phred.
##||MA_peaks_r $dirmacsN_r $dirmacsB_r $qualNum_r


##6c - MAnorm
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
# 3 - diretório de saída para o MACS Narrow
# 4 - diretório de saída para o MACS Broad
# 5 - diretório de saída do MAnorm
# 6 - diretório onde se encontra o programa MAnorm
# 7 - Nome contido no arquivo fq dado a amostra 1 
# 8 - Nome contido no arquivo fq dado a amostra 2
# 9 - Valor de M a ser considerado para as regiões comuns.
# 10 - Valor de M a ser considerado para pico único em cada uma das amostras.
# 11 - Valor máximo de -log10(pvalue) a ser considerado
##|| MAnorm_r $qualNum_r $mergedir $dirmacsN_r $dirmacsB_r $dirMAnorm $dirProgMAnorm $x_tag1_r $x_tag2_r $x_CF_M_Common $x_CF_M_Uniq $x_CF_log10pvalue ;

##7 - IGV-files
##7a - Formatar GFF - FUNÇÃO É IMPORTANTE PARA FORMAT_GFF_R E PARA X_PEAKANALISYS_R.
# 1 - diretório onde se encontra GFF
# 2 - nome do arquivo GFF
##|| format_gff_r $codingdir_r $gffname_r ;

##7b - Busca dos arquivos para IGV
# 1 - diretório para onde serão copiados os arquivos
# 2 - diretório de saída do MAnorm
# 3 - Nome contido no arquivo fq dado a amostra 1 
# 4 - Nome contido no arquivo fq dado a amostra 2
# 5 - diretório de saída para o MACS Narrow
# 6 - diretório de saída para o MACS Broad
# 7 - número a ser considerado para filtro por qualidade phred.
# 8 - diretório Merge
# 9 - diretório onde se encontra GFF
# 10 - nome do arquivo GFF
##|| igv_files_r $dirIGV_r $dirMAnorm $x_tag1_r $x_tag2_r $dirmacsN_r $dirmacsB_r $qualNum_r $mergedir $codingdir_r $gffname_r ;

##8 - Peakanalysis - busca das regiões flanqueadoras e fasta das regiões de controle putativa
# 1 - diretório de saída para os dados
# 2 - diretório onde estão scripts (R_ChiP.R e joinData2.pl
# 3 - diretório onde se encontra GFF
# 4 - nome do arquivo GFF
# 5 - diretório de saída do MAnorm
# 6 - Nome contido no arquivo fq dado a amostra 1 
# 7 - Nome contido no arquivo fq dado a amostra 2
# 8 - diretório do genoma e nome do arquivo.
##|| x_peakanalisys_r $dirAnalisys_r $dirScripts_r $codingdir_r $gffname_r $dirMAnorm $x_tag1_r $x_tag2_r $genome_r;

##9 - Metrics
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
# 3 - diretório de saída pra as métricas


#x_metrics_r $qualNum_r $mergedir $metricsdir $phantomdir $outdirNRFFrip $melngsdir $dirmacsN_r $dirmacsB_r;

#x_SPP_r $qualNum_r $mergedir $phantomdir $metricsdir;

#x_NRF_r $metricsdir $melngsdir;

#x_FRiP_N_r $metricsdir $melngsdir $dirmacsN_r $tag_INPUT_r;

#x_FRiP_B_r $metricsdir $melngsdir $dirmacsB_r $tag_INPUT_r;

##10 - Busca por motivos conservados - Função contém outras funções automatas, fiz assim para facilitar o acréscimo de outras ferramentas de busca por motivos.
motiffinder $qualNum_r $dirmacsN_r $dirmacsB_r $x_tag1_r $x_tag2_r;
