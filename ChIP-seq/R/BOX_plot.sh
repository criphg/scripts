echo $1;
bedtools shuffle -g genome.iff -i $1 > $1"_SHU";
echo $2;
bedtools shuffle -g genome.iff -i $2 > $2"_SHU";

Rscript R_ChIP.R Genes.gff $1 $1;
Rscript R_ChIP.R Genes.gff $2 $2;
Rscript R_ChIP.R Genes.gff $1"_SHU" $1"_SHU";
Rscript R_ChIP.R Genes.gff $2"_SHU" $2"_SHU";

sed "1s/.*/$(head -1 $1".aPeaks" | sed 's/"start"/chr\t"start"/g')/" $1".aPeaks" | awk -F "\t" 'BEGIN{OFS="\t"}{gsub("\"","", $0); print$0}'> temp_X221
mv temp_X221 $1".aPeaks";

sed "1s/.*/$(head -1 $2".aPeaks" | sed 's/"start"/chr\t"start"/g')/" $2".aPeaks" > temp_X221
mv temp_X221 $2".aPeaks";

sed "1s/.*/$(head -1 $1"_SHU.aPeaks" | sed 's/"start"/chr\t"start"/g')/" $1"_SHU.aPeaks" > temp_X221
mv temp_X221 $1"_SHU.aPeaks";

sed "1s/.*/$(head -1 $2"_SHU.aPeaks" | sed 's/"start"/chr\t"start"/g')/" $2"_SHU.aPeaks" > temp_X221
mv temp_X221 $2"_SHU.aPeaks";

Rscript BOXplot_Statistical.R $1".aPeaks" $2".aPeaks" $1"_SHU.aPeaks" $2"_SHU.aPeaks";

