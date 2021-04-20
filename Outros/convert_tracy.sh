

nameout=$(echo $1 | sed "s/.ab1/.fasta/g")
echo "Converting... $1 ...to... $nameout";
var_v=$(echo $1 | sed "s/.ab1//g")
tracy basecall -o temp -f fasta $1
cat temp | awk -v var=$var_v '{if($0 ~ "^>") print">"var ; else print$0}' > $nameout
rm temp;
