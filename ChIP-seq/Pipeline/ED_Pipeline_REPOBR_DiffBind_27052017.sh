####Pipeline corrigido, inclusão do singletons e diretórios de DESEQ e EDGER.


#Pipeline ChIP-seq Rubens 14/07/2015 - revisado 25/05/2017
#Testado - até dia 29/01/2016 - tudo ok até MERGE.

###########IMPORTANTE###########
#NOME DOS ARQUIVOS INPUTS DEVEM SER O MESMO NOME DOS CHIP ATÉ O "_" COM INSERÇÃO DE "-I"  >> ISSO É IMPORTANTE PARA SE DEFINIR NA FUNÇÃO ENRICHMENT O NOME DO CHIP E CONTROLE.
#Ex. se CHIP se chamar SAMPLE1_S1_L001_blabla.fq, o INPUT tem que se chamar SAMPLE1-I_S1_L001_blabla.fq {alternativa, mudar a entrada tagINPUT_r abaixo}
#

#dependências 
#  1 > fastxtoolkit >: sudo apt-get install fastx-toolkit
#  2 > BWA >: sudo apt-get install bwa
#  3 > MACS2 >: https://pypi.python.org/pypi/MACS2/ (more stable version) (version of use macs2 2.1.0.20150731)
#  4 > DiffBind: 
#       Dependencias:
#       sudo apt-get install libcurl4-openssl-dev
#       sudo apt-get install libxml2-dev
#       sudo apt-get install r-cran-xml
#       Entrar no e instalar como sugerido no site http://bioconductor.org/packages/release/bioc/html/DiffBind.html
#       source("https://bioconductor.org/biocLite.R")
#       biocLite("DiffBind")
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
#lista de diretórios contendo fastq.gz para replicata 1
list_r="list_samplesR1";
#lista de diretórios contendo fastq.gz para replicata 2
list_r2="list_samplesR2";


#read length:
read_trim_lgh="75";
#diretórios cotendo fastq.gz dentro da list_samples

##2 - ALIGHMENT
#genome directory and name file
genome_r="Genome/TriTrypDB-9.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta";
#diretório de saída do BAM
outBAM_r="BWA_BAM_REP1";
outBAM_r2="BWA_BAM_REP2";

#número de processos
thread_r="9";

##3 - MERGE
#diretório de saída
mergedir="Merged_REP1";
mergedir2="Merged_REP2";


##4 - Quality filter - phred
#Qualidade phred 
qualNum_r="30";


#########FILTER E MERGER ALLL +++++++ ALTERAR CASO DEIXE DE FAZER ALGUMA DAS ETAPAS DE FILTRO OU REMOÇÃO DE SINGLETOS
mergedirR1=$mergedir"/Q${qualNum_r}/Q${qualNum_r}_pMapped"
mergedirR2=$mergedir2"/Q${qualNum_r}/Q${qualNum_r}_pMapped"

##5 - Enrichment
#diretórios de saída para o MACS
dirmacsN_R1="MACS2_Narrow_Rep1";
dirmacsB_R1="MACS2_Broad_Rep1";
dirmacsN_R2="MACS2_Narrow_Rep2";
dirmacsB_R2="MACS2_Broad_Rep2";


#Tag usada no nome dos arquivos para identificar os inputs
tag_INPUT_r="-I";
#Genome effective size
genomeEffSize_r="17941408";

##6 - DiffBind e ChIPpeakAnno (Junção das duas análises)

#diretório para análise usando DiffBind e ChIPpeakAnno
diffbind_dir="Peak_Analisys";
#diretório onde está os scripts R_ChIP.R e joinData2.pl 
dirScripts_r="scripts";

#dados das amostras para construção da tabela (feito para apenas duas condições, ou dois fatores, ou dois tecidos e uma replicata, para mais de duas amostras ou mais de uma replicata a função deve ser modificada).
df_tissue1="Tcruzi";
df_tissue2="Tcruzi";
df_factorS1="NOTyr";
df_factorS2="NOTyr";
df_condiction1="Resistant";
df_condiction2="Responsive";
df_treatment1="NONE";
df_treatment2="ECM";
df_PeakCallerfile="narrow";
#Tag para arquivos de saída do DiffBind e ChIPpeakAnno
df_output_tag_Narrow="Tcruzi_N";
df_output_tag_Broad="Tcruzi_B";
#Parametro diffbind summits - distancia onde espera-se que estejam os picos a partir do inicio dos fragmentos (IMPORTANTE RODAR SOMENTE O DIFFBIND COM DIFERENTES SUMMITS PARA ESTABELECER PARAMETRO ADEQUADO == MAIOR NUMERO DE REGIOES IDENTIFICADAS)
df_summits="150";  

#Tags de nomeação dos arquivos de entrada DiffBind e IGV-Files (Temporario)
x_tag1_r="Ty";
x_tag2_r="TyM";

##7 - IGV-files
#diretório de saída para os arquivos do IGV
dirIGV_r="IGVfiles";
#diretório onde está o arquivo GFF.
codingdir_r="GFF";
#nome do arquivo GFF.
gffname_r="TriTrypDB-9.0_TcruziCLBrenerEsmeraldo-like_gene.gff";

##8 - Métricas
#diretório de saída pra as métricas
metricsdirR1="metrics_Rep1"; 	#Replicata 1
metricsdirR2="metrics_Rep2"; 	#Replicata 2

#diretório do phantompeakstools
phantomdir="scripts/phantompeakqualtools";
#diretório de saída para cálculo de FRiP e NRF
outdirNRFFrip="NRF_FRiP";
#diretório de localização do mel-ngs
melngsdir="scripts";

##9 - Análise de Motivos
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


##Funções acessórias

STAT_hist_r(){

mkdir STAT;


for i in `ls -d --color=none REP*/*/` ; do

cd $i;
        tag=$(echo $i | awk -F[/-]+ '{print$2}');
        #cout_READS
        echo "Si="$i;
        echo "tag="$tag;

        zcat *R1*fastq.gz | awk 'NR%4==2{print length($0)}' > count_R1_$tag;
        zcat *R2*fastq.gz | awk 'NR%4==2{print length($0)}' > count_R2_$tag;


cd -;

mv $i/count_R* STAT/;


done;


cd STAT;

        cp ../scripts/Hist.R . ;

        #Criação da tabela com tamanho das reads
        for i in `ls --color=none ./count*`; do tag=`echo $i | awk -F"_" '{print$2"_"$3"_"$4}'`; printf "$tag\t$(cat $i | wc -l)\n" ;  done > Read_Size.tab;

        #Criação dos histogramas

        for i in `ls --color=none ./count*`; do Rscript Hist.R $i; done;

cd -;

}



##1 TRIMAGEM

trimag_r () {
echo "##1## TRIMAGEM";
#cat $1;

for a in `cat $1`; do
cd $a;
echo $a;
        for i in `ls --color=auto *fastq.gz`; do
		echo $i;
                zcat $i | fastx_trimmer -Q33 -z -l $2 -i - -o $i"_trim"$2"bp.gz";
        done;
#	pwd;
cd - ;

done;

}

##2 ALIGHMENT
map_r () {
echo "##2## MAPEAMENTO";

if [ ! -f $2".bwt" ] 
	then
	#echo "OK"
	bwa index $2;
fi

#Diretório atual: para voltar no final do loop
local currdir=$(pwd)"/";

for i in `cat $1`;do
cd $i;
        mkdir -p $3;


#Para voltar ao diretório de origem
#        qual=$(echo $i | grep -o "/" | wc -l);
#        dirrec=$(
#        y=$qual;
#        x="-1";
#        while [ $x -lt $y ]
#        do
#        printf "../";
#        x=$[$x+1];
#        done
#        );
#       echo $qual;
#       echo $dirrec;
#       echo ${dirrec}$2;
        ##      


        for a in `ls --color=auto *_R1_*fastq.gz_trim$5bp.gz`;do
           local r2=`echo "$a"|sed 's/_R1_/_R2_/g'`;
           local sam=`echo "$a"|sed 's/_R1_/_R1-2_/g'`;

           bwa aln -q 15 -t $4 -f "./"$3"/"$a"_BWA.sai" ${currdir}$2 $a;
           bwa aln -q 15 -t $4 -f "./"$3"/"$r2"_BWA.sai" ${currdir}$2 $r2;
           bwa sampe -a 750 -f "./"$3"/"$sam"_BWA.sam" ${currdir}$2 "./"$3"/"$a"_BWA.sai" "./"$3"/"$r2"_BWA.sai" $a $r2;
#           samtools view -bSu "./"$3"/"$sam"_BWA.sam"| samtools sort -n -o - samtools_nsort_tmp |samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp |samtools fillmd -u - ${currdir}$2 > "./"$3"/"$sam"_BWA.sam_fixed.bam";

	   #Alteração feita 13/05/2017 para adaptar a nova versão do samtools.
	   samtools view -bSu "./"$3"/"$sam"_BWA.sam"> temp_bam;
           samtools sort -n temp_bam  -o samtools_nsort_tmp; 
	   samtools fixmate samtools_nsort_tmp temp_fix;
	   samtools sort temp_fix -o samtools_csort_tmp;
	   samtools fillmd -u samtools_csort_tmp  ${currdir}$2 > "./"$3"/"$sam"_BWA.sam_fixed.bam";

	   rm temp_bam temp_fix samtools_nsort_tmp samtools_csort_tmp;

           samtools flagstat "./"$3"/"$sam"_BWA.sam_fixed.bam" >> "./"$3"/"$sam"_BWA.sam_fixed.bam_flgst";
	   
        done;
        rm -f *.sam *.sai ;
	

#cd -;
cd ${currdir};

done;


}

map_BWAMEN_r () {
echo "##2## MAPEAMENTO";

local currdir=$(pwd)"/"; 

#Criando Index do Genoma
if [ ! -f $2".bwt" ]
        then
	bwa index $2;
fi


for i in `cat $1`;do
cd $i;
        mkdir -p $3;

        for a in `ls --color=auto *_R1_*fastq.gz_trim$5bp.gz`;do

	echo "##########a = $a"
        echo "##########\$2 = ${dirrec}$2"

           r2=`echo "$a"|sed 's/_R1_/_R2_/g'`;
           sam=`echo "$a"|sed 's/_R1_/_R1-2_/g'`;


           bwa mem -M -t $4 ${currdir}$2 $a $r2 > "./"$3"/"$sam"_BWA.sam";
#           samtools view -bSu "./"$3"/"$sam"_BWA.sam"| samtools sort -n -o - samtools_nsort_tmp |samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp |samtools fillmd -u - ${currdir}$2 > "./"$3"/"$sam"_BWA.sam_fixed.bam";



           #Alteração feita 13/05/2017 para adaptar a nova versão do samtools.
           samtools view -bSu "./"$3"/"$sam"_BWA.sam"> temp_bam;
           samtools sort -n temp_bam  -o samtools_nsort_tmp;
           samtools fixmate samtools_nsort_tmp temp_fix;
           samtools sort temp_fix -o samtools_csort_tmp;
           samtools fillmd -u samtools_csort_tmp  ${currdir}$2 > "./"$3"/"$sam"_BWA.sam_fixed.bam";

           rm temp_bam temp_fix samtools_nsort_tmp.bam samtools_csort_tmp.bam;



           samtools flagstat "./"$3"/"$sam"_BWA.sam_fixed.bam" >> "./"$3"/"$sam"_BWA.sam_fixed.bam_flgst";
        done;
        rm -f *.sam *.sai ;
cd ${currdir};
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
        samtools sort $Name -o $1"/"$OutName".bam";
        samtools index $1"/"$OutName".bam";
	samtools flagstat $1"/"$OutName".bam" > $1"/"$OutName"_flgst";

done;

rm -f $1"/*fixed.bam";

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
                samtools sort $Name -o $Name2".bam";
                samtools index $Name2".bam";
                samtools flagstat $Name2".bam" > $Name2"_flgst";
        done;

        rm -f Q$1/*_tempX;

cd - ;

}


##4.b QUALITY FILTER PROPER PAIR AND BOTH PAIRS UNIQUE MAPPED.
################FALTA EDITAR!!!!! COLOCAR JUNTO COM O DE CIMA E EDITAR ENRICHMENT
#Acreditei no flagstat... possível problema com a seleção dos pares... verificando.
# Filtro dos arquivos de  qualidade Q30, por propriamente mapeados -f 3, onde ambos os membros do par são unicamente mapeados.

filter_mateBWA_XTAR_r() {

echo "###4.b### QUALITY PAIR";

cd $2"/Q"$1 ;

mkdir "Q"$1"_pMapped";

        for i in `ls --color=auto Q${1}_*index.bam`; do

        cd "Q"$1"_pMapped";
                outname=`echo $i | awk -F"." '{print$1}'`;
                echo "Outname =" $outname;
                samtools view -H "../"$i > temp_file0;  #BAM HEADER
                samtools view -f 3 "../"$i > "temp_file" ;
#                cat temp_file | awk '{print $1}' | sort | uniq -u > list_temp;                    #Proper pair verdadeiros
                cat temp_file | grep -v "^@" | awk '{if ($12=="XT:A:R") print$1}' > list_temp;   #Proper pair selecionado pela tag

		cat temp_file | fgrep -w -v -f list_temp > temp_file00; #Remoção do list_temp
                cat temp_file0 temp_file00 > temp_file2;		#Adição de cabeçalho


		samtools view -bS temp_file2 > temp_file3
                samtools sort temp_file3 -o $outname"_pMapped_index.bam";
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

for i in `ls --color=auto ${3}/Q${2}_*index.bam | grep '\\'${4}`; do
	echo $i;
        control=`echo "../$i"`;
        chip=`echo "../"$i | sed 's/\-I//g'`;
        name=`echo $chip | awk -F"/" '{print$NF}' | awk -F"." '{print$1}'`;

#name=`echo $name | awk -F"_" '{print$2}'`;

        echo "Control==="$control;
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

for i in `ls --color=auto ${3}/Q${2}_*index.bam | grep '\\'${4}`; do
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




##6 - DiffBind <<<--
##6.1 - Criação da tabela de dados sobre amostras (csv).






#função complementar ao DiffBind, ajustando os cabeçalhos das tabelas



Head_Form_R(){

        ###Formatando arquivos de saida DESEQ2
cd $1;

        echo "DiffBind head output file formating";
        ##Arrumando o cabeçalho do arquivo completo
        #Narrow
        awk 'NR==1 {OFS = "\t"; $1 = "ID\t" $1; print $0} NR!=1 {print$0}' ${df_output_tag_Narrow}_DB_analisys.out > ${df_output_tag_Narrow}_DB_analisys.peaks
        #Broad
        awk 'NR==1 {OFS = "\t"; $1 = "ID\t" $1; print $0} NR!=1 {print$0}' ${df_output_tag_Broad}_DB_analisys.out > ${df_output_tag_Broad}_DB_analisys.peaks

        ##Separando picos de Resistent e Responsive (MIN_FC = 2, FDR < 0.05)
        #Narrow
        awk 'NR==1 {print $0} NR!=1 {if ($10 > 0) print $0}' ${df_output_tag_Narrow}_DB_analisys.peaks > ${df_output_tag_Narrow}_Resistent.peaks;
	awk 'NR==1 {print $0} NR!=1 {if ($10 < 0) print $0}' ${df_output_tag_Narrow}_DB_analisys.peaks > ${df_output_tag_Narrow}_Responsive.peaks;


        #Broad
        awk 'NR==1 {print $0} NR!=1 {if ($10 > 0) print $0}' ${df_output_tag_Broad}_DB_analisys.peaks > ${df_output_tag_Broad}_Resistent.peaks;
        awk 'NR==1 {print $0} NR!=1 {if ($10 < 0) print $0}' ${df_output_tag_Broad}_DB_analisys.peaks > ${df_output_tag_Broad}_Responsive.peaks;

        #Removendo os arquivos de saida não formatados
 #       rm -f *.out;

cd - ;

echo "ChIPpeakAnno head output file formating";
cd $2;

        ##Arrumando cabeçalho do arquivo de saído do ChIPpeakAnno
        #Narrow
        awk 'NR==1 {OFS = "\t"; $1 = "ID\t" $1; print $0} NR!=1 {print$0}' ${df_output_tag_Narrow}.out > temp_N.out2;
        #Acrescentando informacões com o nome dos genes no Tritryp e ID
        perl ../../${dirScripts_r}/joinData2.pl ../../${codingdir_r}/Genes.gff temp_N.out2 > ${df_output_tag_Narrow}.aPeaks;
        
        #Broad
        awk 'NR==1 {OFS = "\t"; $1 = "ID\t" $1; print $0} NR!=1 {print$0}' ${df_output_tag_Broad}.out > temp_B.out2;
        #Acrescentando informacões com o nome dos genes no Tritryp e ID
        perl ../../${dirScripts_r}/joinData2.pl ../../${codingdir_r}/Genes.gff temp_B.out2 > ${df_output_tag_Broad}.aPeaks;

        #Removendo os arquivos de saida não formatados
        rm -f *.out2;


cd - ;
        

}


diffbindanalisys(){
echo "Criando tabela para análise DiffBind...";


Begin_D=`pwd`;
echo $Begin_D;

mkdir -p $1;
cd $1;
	#Tabela Narrow:
	#Cabeçalho:
	echo "SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller" > tabdiff_N.csv;
	echo "SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller" > tabdiff_B.csv;


	#Amostra1 - REP1
	BAMREP1CHIP=`ls --color=auto ../${mergedir}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag1_r}_*.bam`;
	BAMREP1CTRL=`ls --color=auto ../${mergedir}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag1_r}${tag_INPUT_r}_*.bam`;
	BEDREP1_N=`ls --color=auto ../$dirmacsN_R1/Q${qualNum_r}*${x_tag1_r}_*.narrowPeak`;
	BEDREP1_B=`ls --color=auto ../$dirmacsB_R1/Q${qualNum_r}*${x_tag1_r}_*.broadPeak`;
	echo "BAM1==$BAMREP1CHIP";
	echo "BAM2==$BAMREP1CTRL";
	echo "BED_N==$BEDREP1_N";
	echo "BED_B==$BEDREP1_B";	

#tabela Narrow
	echo "${x_tag1_r}1,$df_tissue1,$df_factorS1,$df_condiction1,$df_treatment1,1,$BAMREP1CHIP,${x_tag1_r}${tag_INPUT_r}1,$BAMREP1CTRL,$BEDREP1_N,$df_PeakCallerfile" >> tabdiff_N.csv;
#tabela Broad
        echo "${x_tag1_r}1,$df_tissue1,$df_factorS1,$df_condiction1,$df_treatment1,1,$BAMREP1CHIP,${x_tag1_r}${tag_INPUT_r}1,$BAMREP1CTRL,$BEDREP1_B,$df_PeakCallerfile" >> tabdiff_B.csv;
	

	#Amostra1 - REP2
        BAMREP2CHIP=`ls --color=auto ../${mergedir2}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag1_r}_*.bam`;
        BAMREP2CTRL=`ls --color=auto ../${mergedir2}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag1_r}${tag_INPUT_r}_*.bam`;
        BEDREP2_N=`ls --color=auto ../$dirmacsN_R2/Q${qualNum_r}*${x_tag1_r}_*.narrowPeak`;
        BEDREP2_B=`ls --color=auto ../$dirmacsB_R2/Q${qualNum_r}*${x_tag1_r}_*.broadPeak`;
        echo "BAM1==$BAMREP2CHIP";
        echo "BAM2==$BAMREP2CTRL";      
        echo "BED_N==$BEDREP2_N";
        echo "BED_B==$BEDREP2_B";
    
#tabela Narrow
        echo "${x_tag1_r}2,$df_tissue1,$df_factorS1,$df_condiction1,$df_treatment1,2,$BAMREP2CHIP,${x_tag1_r}${tag_INPUT_r}2,$BAMREP2CTRL,$BEDREP2_N,$df_PeakCallerfile" >> tabdiff_N.csv
#tabela Broad
        echo "${x_tag1_r}2,$df_tissue1,$df_factorS1,$df_condiction1,$df_treatment1,2,$BAMREP2CHIP,${x_tag1_r}${tag_INPUT_r}2,$BAMREP2CTRL,$BEDREP2_B,$df_PeakCallerfile" >> tabdiff_B.csv


	#Amostra2 - REP1
	BAMREP1CHIP=`ls --color=auto ../${mergedir}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag2_r}_*.bam`;
        BAMREP1CTRL=`ls --color=auto ../${mergedir}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag2_r}${tag_INPUT_r}_*.bam`;
        BEDREP1_N=`ls --color=auto ../$dirmacsN_R1/Q${qualNum_r}*${x_tag2_r}_*.narrowPeak`;
        BEDREP1_B=`ls --color=auto ../$dirmacsB_R1/Q${qualNum_r}*${x_tag2_r}_*.broadPeak`;
        echo "BAM1==$BAMREP1CHIP";
        echo "BAM2==$BAMREP1CTRL";
        echo "BED_N==$BEDREP1_N";
	echo "BED_B==$BEDREP1_B";

#tabela Narrow
        echo "${x_tag2_r}1,$df_tissue2,$df_factorS2,$df_condiction2,$df_treatment2,1,$BAMREP1CHIP,${x_tag2_r}${tag_INPUT_r}1,$BAMREP1CTRL,$BEDREP1_N,$df_PeakCallerfile" >> tabdiff_N.csv;
#tabela Broad
        echo "${x_tag2_r}1,$df_tissue2,$df_factorS2,$df_condiction2,$df_treatment2,1,$BAMREP1CHIP,${x_tag2_r}${tag_INPUT_r}1,$BAMREP1CTRL,$BEDREP1_B,$df_PeakCallerfile" >> tabdiff_B.csv;


	#Amostra2 - REP2
	BAMREP2CHIP=`ls --color=auto ../${mergedir2}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag2_r}_*.bam`;
        BAMREP2CTRL=`ls --color=auto ../${mergedir2}/Q${qualNum_r}/Q${qualNum_r}_pMapped/Q${qualNum_r}*${x_tag2_r}${tag_INPUT_r}_*.bam`;
        BEDREP2_N=`ls --color=auto ../$dirmacsN_R2/Q${qualNum_r}*${x_tag2_r}_*.narrowPeak`;
        BEDREP2_B=`ls --color=auto ../$dirmacsB_R2/Q${qualNum_r}*${x_tag2_r}_*.broadPeak`;
        echo "BAM1==$BAMREP2CHIP";
        echo "BAM2==$BAMREP2CTRL";
        echo "BED_N==$BEDREP2_N";
        echo "BED_B==$BEDREP2_B";

#tabela Narrow
        echo "${x_tag2_r}2,$df_tissue2,$df_factorS2,$df_condiction2,$df_treatment2,2,$BAMREP2CHIP,${x_tag2_r}${tag_INPUT_r}2,$BAMREP2CTRL,$BEDREP2_N,$df_PeakCallerfile" >> tabdiff_N.csv
#tabela Broad
        echo "${x_tag2_r}2,$df_tissue2,$df_factorS2,$df_condiction2,$df_treatment2,2,$BAMREP2CHIP,${x_tag2_r}${tag_INPUT_r}2,$BAMREP2CTRL,$BEDREP2_B,$df_PeakCallerfile" >> tabdiff_B.csv

#Iniciando DiffBind e ChIPpeakAnno
	
	echo "Starting R script, DiffBind e ChIPpeakAnno";
        Rscript ../${dirScripts_r}/DiffBind_ChIPpeakAnno_DESEQ2.R ../${codingdir_r}/Genes.gff tabdiff_N.csv ${df_output_tag_Narrow} $df_summits;
        Rscript ../${dirScripts_r}/DiffBind_ChIPpeakAnno_DESEQ2.R ../${codingdir_r}/Genes.gff tabdiff_B.csv ${df_output_tag_Broad} $df_summits;
        Rscript ../${dirScripts_r}/DiffBind_ChIPpeakAnno_EDGER.R ../${codingdir_r}/Genes.gff tabdiff_N.csv ${df_output_tag_Narrow} $df_summits;
        Rscript ../${dirScripts_r}/DiffBind_ChIPpeakAnno_EDGER.R ../${codingdir_r}/Genes.gff tabdiff_B.csv ${df_output_tag_Broad} $df_summits;

        #Uso da função para ajuste de cabeçalho # postada anterior a essa.
        Head_Form_R DiffBind_D_DESEQ2 ChIPpeakAnno_D_DESEQ2;
        Head_Form_R DiffBind_D_EDGER ChIPpeakAnno_D_EDGER;


	
cd .. ;


}


##ARQUIVO GFF PRECISA SER PREPARADO. ELIMINAÇÃO DO CABEÇALHO E SELEÇÃO SOMENTE DOS GENES.

format_gff_r(){
cd ${codingdir_r};

	grep -v "^#" ${gffname_r} | awk '{if($3 == "gene") print$0}' > Genes.gff;

cd .. ;
}



igv_files_r(){

############## IGV-files

echo "###7### IGV-files";


echo "##Copiando e renomeando arquivos para o diretório ${dirIGV_r} ##";

mkdir ${dirIGV_r};
curr_dir=`pwd .`;

cd ${dirIGV_r};


#########EDITAR ARQUIVOS PARA ADMIIR QUALQUER TAG

        mkdir Broad;
        mkdir Narrow;


	#Diffbind
	## convertendo peaks para bed e wig
	
	#Narrow
	#all_bed and bedgraph
	cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Narrow}_DB_analisys.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4}' > Narrow/${df_output_tag_Narrow}.bed;
	cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Narrow}_DB_analisys.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4"\t"$10}' > Narrow/${df_output_tag_Narrow}.bedGraph;

	#Uniq_Ty
	cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Narrow}_Resistent.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4}' > Narrow/${x_tag1_r}_Uniq_N.bed;
	cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Narrow}_Resistent.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4"\t"$10}' > Narrow/${x_tag1_r}_Uniq_N.bedGraph;

	#Uniq_TyM
	cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Narrow}_Responsive.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4}' > Narrow/${x_tag2_r}_Uniq_N.bed;
        cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Narrow}_Responsive.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4"\t"$10}' > Narrow/${x_tag2_r}_Uniq_N.bedGraph;




#narrow_Rep1
        cp $curr_dir/${dirmacsN_R1}/{*.bdg,*.narrowPeak} ./Narrow;
        mv ./Narrow/Q${qualNum_r}_${x_tag1_r}_MERGED*_pileup.bdg ./Narrow/ChIP_${x_tag1_r}_MACS_N_Rep1.bdg;
        mv ./Narrow/Q${qualNum_r}_${x_tag2_r}_MERGED*_pileup.bdg ./Narrow/ChIP_${x_tag2_r}_MACS_N_Rep1.bdg;
        mv ./Narrow/Q${qualNum_r}_${x_tag1_r}_MERGED*_lambda.bdg ./Narrow/INPUT_${x_tag1_r}_MACS_N_Rep1.bdg;
        mv ./Narrow/Q${qualNum_r}_${x_tag2_r}_MERGED*_lambda.bdg ./Narrow/INPUT_${x_tag2_r}_MACS_N_Rep1.bdg;

        mv ./Narrow/Q${qualNum_r}_${x_tag1_r}_MERGED*.narrowPeak ./Narrow/${x_tag1_r}_MACS_N_Rep1.narrowPeak;
        mv ./Narrow/Q${qualNum_r}_${x_tag2_r}_MERGED*.narrowPeak ./Narrow/${x_tag2_r}_MACS_N_Rep1.narrowPeak;

#narrow_Rep2
        cp $curr_dir/${dirmacsN_R2}/{*.bdg,*.narrowPeak} ./Narrow;
        mv ./Narrow/Q${qualNum_r}_${x_tag1_r}_MERGED*_pileup.bdg ./Narrow/ChIP_${x_tag1_r}_MACS_N_Rep2.bdg;
        mv ./Narrow/Q${qualNum_r}_${x_tag2_r}_MERGED*_pileup.bdg ./Narrow/ChIP_${x_tag2_r}_MACS_N_Rep2.bdg;
        mv ./Narrow/Q${qualNum_r}_${x_tag1_r}_MERGED*_lambda.bdg ./Narrow/INPUT_${x_tag1_r}_MACS_N_Rep2.bdg;
        mv ./Narrow/Q${qualNum_r}_${x_tag2_r}_MERGED*_lambda.bdg ./Narrow/INPUT_${x_tag2_r}_MACS_N_Rep2.bdg;

	mv ./Narrow/Q${qualNum_r}_${x_tag1_r}_MERGED*.narrowPeak ./Narrow/${x_tag1_r}_MACS_N_Rep2.narrowPeak;
        mv ./Narrow/Q${qualNum_r}_${x_tag2_r}_MERGED*.narrowPeak ./Narrow/${x_tag2_r}_MACS_N_Rep2.narrowPeak;

        cp $curr_dir/${codingdir_r}/Genes.gff ./Narrow;



	#Broad
	#all_bed and bedgraph
        cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Broad}_DB_analisys.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4}' > Broad/${df_output_tag_Broad}.bed;
        cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Broad}_DB_analisys.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4"\t"$10}' > Broad/${df_output_tag_Broad}.bedGraph;

        #Uniq_Ty
        cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Broad}_Resistent.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4}' > Broad/${x_tag1_r}_Uniq_B.bed;
        cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Broad}_Resistent.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4"\t"$10}' > Broad/${x_tag1_r}_Uniq_B.bedGraph;

        #Uniq_TyM
        cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Broad}_Responsive.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4}' > Broad/${x_tag2_r}_Uniq_B.bed;
        cat $curr_dir/${diffbind_dir}/DiffBind_D_DESEQ2/${df_output_tag_Broad}_Responsive.peaks | tail -n +2 | sed 's/\"//g' | awk '{print $2"\t"$3"\t"$4"\t"$10}' > Broad/${x_tag2_r}_Uniq_B.bedGraph;


#Broad_Rep1
        cp $curr_dir/${dirmacsB_R1}/{*.bdg,*.broadPeak} ./Broad;
        mv ./Broad/Q${qualNum_r}_${x_tag1_r}_MERGED*_pileup.bdg ./Broad/ChIP_${x_tag1_r}_MACS_B_Rep1.bdg;
        mv ./Broad/Q${qualNum_r}_${x_tag2_r}_MERGED*_pileup.bdg ./Broad/ChIP_${x_tag2_r}_MACS_B_Rep1.bdg;
        mv ./Broad/Q${qualNum_r}_${x_tag1_r}_MERGED*_lambda.bdg ./Broad/INPUT_${x_tag1_r}_MACS_B_Rep1.bdg;
        mv ./Broad/Q${qualNum_r}_${x_tag2_r}_MERGED*_lambda.bdg ./Broad/INPUT_${x_tag2_r}_MACS_B_Rep1.bdg;

	mv ./Broad/Q${qualNum_r}_${x_tag1_r}_MERGED*.broadPeak ./Broad/${x_tag1_r}_MACS_B_Rep1.broadPeak;
        mv ./Broad/Q${qualNum_r}_${x_tag2_r}_MERGED*.broadPeak ./Broad/${x_tag2_r}_MACS_B_Rep1.broadPeak;

#Broad_Rep2
	cp $curr_dir/${dirmacsB_R2}/{*.bdg,*.broadPeak} ./Broad;
        mv ./Broad/Q${qualNum_r}_${x_tag1_r}_MERGED*_pileup.bdg ./Broad/ChIP_${x_tag1_r}_MACS_B_Rep2.bdg;
        mv ./Broad/Q${qualNum_r}_${x_tag2_r}_MERGED*_pileup.bdg ./Broad/ChIP_${x_tag2_r}_MACS_B_Rep2.bdg;
        mv ./Broad/Q${qualNum_r}_${x_tag1_r}_MERGED*_lambda.bdg ./Broad/INPUT_${x_tag1_r}_MACS_B_Rep2.bdg;
        mv ./Broad/Q${qualNum_r}_${x_tag2_r}_MERGED*_lambda.bdg ./Broad/INPUT_${x_tag2_r}_MACS_B_Rep2.bdg;

	mv ./Broad/Q${qualNum_r}_${x_tag1_r}_MERGED*.broadPeak ./Broad/${x_tag1_r}_MACS_B_Rep2.broadPeak;
        mv ./Broad/Q${qualNum_r}_${x_tag2_r}_MERGED*.broadPeak ./Broad/${x_tag2_r}_MACS_B_Rep2.broadPeak;

	cp $curr_dir/${codingdir_r}/Genes.gff ./Broad;

cd $curr_dir ;
}


##9 - Métricas


#Funções acessórias para métricas vindas do MAnorm.

# Tranformar BAMtoBED dos bam pós filtros para SPP.

BAMtoBED_MA_reads_r() {
echo "##Getting Reads BAMtoBED##"
echo"";

if [ "${4}" == "XXPROCEED" ]; then
	cd ${2}
else
	cd "${2}/Q${1}/Q${1}_pMapped/"
fi

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

if [ ! "$4" == "XXPROCEED" ]; then
cd ../../../;
else
cd ../;
fi

}

BAMtoBED_MA_reads_INPUTS_r() {
echo "##Getting Reads BAMtoBED_INPUT##"
echo"";

if [ "$4" == "XXPROCEED" ]; then
cd ${2};
else
cd "${2}/Q${1}/Q${1}_pMapped/";
fi

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
	
if [ ! "$4" == "XXPROCEED" ]; then
cd ../../../;
else
cd ../;
fi
}




###########################
#########################Função desnecessária
###############################

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



#Função SPP
#formatar arquivo SPP localizado na pasta /media/sf_D_DRIVE/Dados/Tcruzi/SPP_Analise , arquivo SPP_quality_analisys.R
#verificar uso do runSPP.R para calculo dessa métrica /Tools/phantompeakqualtools
#Fazer o link simbólico do Ty.bed para Ty.tagAlign
#verificar outras métricas em ~/Tools/mel-ngs/mel-chipseq/chipseq-metrics 


#Não é necessário separar para broad e narrow, SPP leve em consideração somente reads mapeadas para buscar Cross-correlation.
x_SPP_r () {

#cria diretório de métricas, se já existir não mostra erro.
mkdir -p ${4};
curr_dir=`pwd .`;

cd ${4};


	echo "##9## Metrics - SPP - Cross-correlation";

	#Busca dos arquivos bed, link simbólico mudando o nome do arquivo para .tagAlign >> necessário para o SPP
	for i in `ls $curr_dir/${2}/Q${1}/Q${1}_pMapped/{bed,bed/INPUT}/read*.bed`; do
        	namef=`echo $i | awk -F"/" '{print$NF}' | awk -F"." '{print$1}' | awk -F"reads_" '{print$2}'`;
        	echo $namef;
        	ln -s $i ./${namef}".tagAlign" ;
	done;


	mkdir SPP_CrossCorrelation;
	cd SPP_CrossCorrelation;

        	#Busca SPP.R

        	ln -s $curr_dir/${3}/run_spp.R .;


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



##10 - Apresentação dos resultados
##10.1 - Genomation

x_genomation_r(){

mkdir -p Heatmap_genomation;

echo "###10### - Results presentation";
echo "###10.1### - Heatmap genomation";

cd Heatmap_genomation;

	Rscript ../scripts/genomation_heatmap.R $diffbind_dir"/DiffBind_D_DESEQ2/" $mergedir"/Q"$qualNum_r"/Q"$qualNum_r"_pMapped" $mergedir2"/Q"$qualNum_r"/Q"$qualNum_r"_pMapped/" $x_tag1_r $x_tag2_r ;

cd .. ;

}

##11 - Busca de motivos 



##11a - função HOMER



x_homer_r () {

echo "Motif searching - HOMER";
mkdir -p $1;
cd $1;	
	#transformar o arquivo do MAnorm S1_UNIQ ou S2_UNIQ em bed
	cat ../../$2/$3/S1_UNIQ_PEAK.out | awk -F"\t" 'BEGIN{OFS="\t"} {print$2,$3,$4,$1}' > S1_UNIQ_PEAK_${3}.bed;
        cat ../../$2/$3/S2_UNIQ_PEAK.out | awk -F"\t" 'BEGIN{OFS="\t"} {print$2,$3,$4,$1}' > S2_UNIQ_PEAK_${3}.bed;

	#limha de comando para homer
	#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
	findMotifsGenome.pl S1_UNIQ_PEAK_${3}.bed "../../"$7 $3'/'$6 -size $4
	findMotifsGenome.pl S2_UNIQ_PEAK_${3}.bed "../../"$7 $3'/'$6 -size $5

cd ..;

}


motiffinder () {
echo "###11### - Motif Searching";

mkdir -p $motifdir;
cd $motifdir;
        #busca por domínio com parametros definidos para Broad e Narrow.
        x_homer_r $homerdir $dirMAnorm Narrow $homerfragsize_S1 $homerfragsize_S2 $homeroutdir $genome_r;
        x_homer_r $homerdir $dirMAnorm Broad $homerfragsize_S1 $homerfragsize_S2 $homeroutdir $genome_r;



	#busca por domínio com size médio definido pelo MACS para Broad e Narrow.
	#size médio MACS

	fragBroad_Sample1=`cat ../${3}/Q${1}_${4}_*_index_b.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
	fragBroad_Sample2=`cat ../${3}/Q${1}_${5}_*_index_b.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
	fragNarrow_Sample1=`cat ../${2}/Q${1}_${4}_*_index_n.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;
	fragNarrow_Sample2=`cat ../${2}/Q${1}_${5}_*_index_n.log | grep "fragment size =" | awk -F"=" '{print$2}' | tr -d ' '` ;


#$qualNum_r $dirmacsN_r $dirmacsB_r $x_tag1_r $x_tag2_r


	x_homer_r $homerdir $dirMAnorm Narrow $fragNarrow_Sample1 $fragNarrow_Sample2 Result_HOMER_FSizeMACS $genome_r;
        x_homer_r $homerdir $dirMAnorm Broad $fragBroad_Sample1 $fragBroad_Sample2 Result_HOMER_FSizeMACS $genome_r;




cd .. ;

}


#=============================================================================
#CHAMADA DE FUNÇÕES:

#####Statistica e Histogramas
##STAT_hist_r ;

##1 - TRIMAGEM

#deve ser informado: 
# 1 - list_samples com os diretórios de cada amostra contendo os fastq.gz!
# 2 - tamanho da região a ser deixada.
#trimag_r $list_r $read_trim_lgh; trimag_r $list_r2 $read_trim_lgh;

##2 - ALIGHMENT

#deve ser informado:
# 1 - list_samples com os diretórios de cada amostra contendo os fastq.gz!
# 2 - diretório do genoma e nome do arquivo.
# 3 - diretório de saída do bam

#BWA_SW
#map_r $list_r $genome_r $outBAM_r $thread_r $read_trim_lgh; 
#map_r $list_r2 $genome_r $outBAM_r2 $thread_r $read_trim_lgh;

#BWA_MEN
#map_BWAMEN_r $list_r $genome_r $outBAM_r $thread_r $read_trim_lgh;
#map_BWAMEN_r $list_r2 $genome_r $outBAM_r2 $thread_r $read_trim_lgh;

##3 - MERGE
#deve ser informado:
# 1 - ditetório de saida para o Merge
#merge_r $mergedir $list_r $outBAM_r;
#merge_r $mergedir2 $list_r2 $outBAM_r2;


##4 - Quality Filer Phred
#deve ser informado:
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
#quality_pherd_r $qualNum_r $mergedir;
#quality_pherd_r $qualNum_r $mergedir2;

##4b - Filter mates with r1 XT:A:R, r2 XT:A:U
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
#filter_mateBWA_XTAR_r $qualNum_r $mergedir;
#filter_mateBWA_XTAR_r $qualNum_r $mergedir2;

##5 - Enrichment
##NARROW
# 1 - diretório de saída para o MACS
# 2 - número a ser considerado para filtro por qualidade phred.
# 3 - diretório Merge
# 4 - Tag usada nos arquivos de INPUT para o MACS diferenciar INPUTS de CHIP.
# 5 - Genome effective size 
#x_macs2Narrow_r $dirmacsN_R1 $qualNum_r $mergedirR1 $tag_INPUT_r $genomeEffSize_r ;
#x_macs2Narrow_r $dirmacsN_R2 $qualNum_r $mergedirR2 $tag_INPUT_r $genomeEffSize_r ;

##BROAD
# 1 - diretório de saída para o MACS
# 2 - número a ser considerado para filtro por qualidade phred.
# 3 - diretório Merge
# 4 - Tag usada nos arquivos de INPUT para o MACS diferenciar INPUTS de CHIP.
# 5 - Genome effective size 
#x_macs2Broad_r $dirmacsB_R1 $qualNum_r $mergedirR1 $tag_INPUT_r $genomeEffSize_r ;
#x_macs2Broad_r $dirmacsB_R2 $qualNum_r $mergedirR2 $tag_INPUT_r $genomeEffSize_r ;

##6 - DiffBind
# 1 - diretório de saida para o DiffBind
diffbindanalisys $diffbind_dir;

##7 - IGV-files
##7a - Formatar GFF - FUNÇÃO É IMPORTANTE PARA FORMAT_GFF_R E PARA X_PEAKANALISYS_R.
# 1 - diretório onde se encontra GFF
# 2 - nome do arquivo GFF
#format_gff_r ;

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
#igv_files_r ;

##8 - Metrics - Adequação do arquivo contendo reads de Ty e TyM para calculo de Métricas
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
#BAMtoBED_MA_reads_r $qualNum_r $mergedir $tag_INPUT_r;
#BAMtoBED_MA_reads_r $qualNum_r $mergedir2 $tag_INPUT_r;

##6a.2 - Obtenção dos beds também para os INPUTS. Necessário para cálculos de métricas ou mesmo para vizualização no IGV.
#BAMtoBED_MA_reads_INPUTS_r $qualNum_r $mergedir $tag_INPUT_r;
#BAMtoBED_MA_reads_INPUTS_r $qualNum_r $mergedir2 $tag_INPUT_r;

##9.1 - Metrics
# 1 - número a ser considerado para filtro por qualidade phred.
# 2 - diretório Merge
# 3 - diretório de saída pra as métricas


#x_metrics_r $qualNum_r $mergedir $metricsdir $phantomdir $outdirNRFFrip $melngsdir $dirmacsN_r $dirmacsB_r;

#x_SPP_r $qualNum_r $mergedir $phantomdir $metricsdirR1;
#x_SPP_r $qualNum_r $mergedir2 $phantomdir $metricsdirR2;

#x_NRF_r $metricsdirR1 $melngsdir;
#x_NRF_r $metricsdirR2 $melngsdir;

#x_FRiP_N_r $metricsdirR1 $melngsdir $dirmacsN_R1 $tag_INPUT_r;
#x_FRiP_N_r $metricsdirR2 $melngsdir $dirmacsN_R2 $tag_INPUT_r;

#x_FRiP_B_r $metricsdirR1 $melngsdir $dirmacsB_R1 $tag_INPUT_r;
#x_FRiP_B_r $metricsdirR2 $melngsdir $dirmacsB_R2 $tag_INPUT_r;


##10 - Apresentação dos dados
##10.1 - Heatmap de enriquecimento - genomation

#x_genomation_r 






##11 - Busca por motivos conservados - Função contém outras funções automatas, fiz assim para facilitar o acréscimo de outras ferramentas de busca por motivos.
#||||motiffinder $qualNum_r $dirmacsN_r $dirmacsB_r $x_tag1_r $x_tag2_r;
