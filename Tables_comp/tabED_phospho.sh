#tabED.sh <INPUT-TAB> <OUTPUT>


echo -e "Pept\tGeneID\tOrganism\tPtn_Name\tChrm_Location\tlenght\tSequence_SO\tGO\tAverage_Ty\tAverage_TyM\tRation_TyM/Ty\tp_Value\tlog_ration\tlog_neg_pvalue" > $2;
tail -n +2 $1 >> $2; 

sed -i -e "s/ | /\t/g" $2;

sed -i  "s/,/\t/g" $2;

sed -i  "s/organism=//g" $2;

sed -i -e "s/product=//g" $2;

sed -i -e "s/sequence_SO=//g" $2 ;

sed -i -e "s/SO=//g" $2 ;

sed -i -e "s/location=//g" $2 


sed -i -e "s/length=//g" $2 
