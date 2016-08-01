#Edited by criph 10/03/2016

#!/bin/bash
if [ $# -ne 3 ]
then
  echo "Usage: `basename $0` fold-change_cutoff_unbiased fold-change_cutoff_biased -log10p_cutoff_biased"
  exit
fi

#editado por mim para acrescer o valor de M no arquivo final gerado, acrescimo de $7.
#sed '1d' MAnorm_result.xls | awk 'BEGIN {OFS="\t"}{print $1,$2,$3>"MAnorm_result.tmp"}'
sed '1d' MAnorm_result.xls | awk 'BEGIN {OFS="\t"}{$7=sprintf("%.14f", $7); print$0}' | sort -k7,7nr -k8,8nr > MAnorm_result.tmp;

#awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($7 >= var1 && $9 > var2)print $1,$2,$3,$7>"sample1_uniq_peaks.bed"}'  MAnorm_result.xls
#awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($7 >= var1 && $9 > var2) {sum+=1; print $1,$2,$3,"S1_Peak"sum,$7>"sample1_uniq_peaks.bed"}}'  MAnorm_result.tmp;

(head -1 MAnorm_result.xls | sed 's/chr/PeakID\tchr/g'; awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($7 >= var1 && $9 > var2) {sum+=1; $1="S1_Peak"sum"\t"$1; print$0}}'  MAnorm_result.tmp) > "S1_UNIQ_PEAK.out"
sed '1d' S1_UNIQ_PEAK.out | awk 'BEGIN{OFS = "\t"}{print $2,$3,$4,$1,$8}' > "sample1_uniq_peaks.bed"


cat MAnorm_result.tmp | sort -k7,7n -k8,8n > MAnorm_result.tmp2;


#awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($7 <= -var1 && $9 > var2)print $1,$2,$3,$7>"sample2_uniq_peaks.bed"}'  MAnorm_result.xls
#awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($7 <= -var1 && $9 > var2){sum+=1; print $1,$2,$3,"S2_Peak"sum,$7>"sample2_uniq_peaks.bed"}}'  MAnorm_result.tmp2;

(head -1 MAnorm_result.xls | sed 's/chr/PeakID\tchr/g'; awk -v var1=$2 -v var2=$3 'BEGIN {OFS="\t"} {if($7 <= -var1 && $9 > var2) {sum+=1; $1="S2_Peak"sum"\t"$1; print$0}}'  MAnorm_result.tmp2) > "S2_UNIQ_PEAK.out"
sed '1d' S2_UNIQ_PEAK.out | awk 'BEGIN{OFS = "\t"}{print $2,$3,$4,$1,$8}' > "sample2_uniq_peaks.bed"



#awk -v var=$1 'BEGIN {OFS="\t"} {if($7<var && $7>-var)print $1,$2,$3,$7>"unbiased_peaks_ED.tmp"}'  MAnorm_result.xls
awk -v var=$1 'BEGIN {OFS="\t"} {if($7<var && $7>-var){sum+=1; print $1,$2,$3,"Common_Peak"sum,$7>"unbiased_peaks_ED.tmp"}}'  MAnorm_result.tmp;

#editado por mim, para evitar erro no mergeBed, também modificado arquivo unbiased_peaks.tmp para unbiased_peaks_ED.tmp linha 12. Adicionado rm unbiased_peaks_ED.tmp linha 21.

cat unbiased_peaks_ED.tmp | sort -k1,1 -k2,2n -k3,3n > unbiased_peaks.tmp;

mergeBed -i unbiased_peaks.tmp > unbiased_peaks.bed

rm unbiased_peaks.tmp
rm unbiased_peaks_ED.tmp
rm MAnorm_result.tmp
rm MAnorm_result.tmp2
