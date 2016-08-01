#baixa as sequencias fasta do ensembl a partir de uma lista contendo ensembl ID.


for i in `cat Gene_ensID_Tabela_completa_Edger`; do echo -e ">"$i; echo `elinks --dump http://rest.ensembl.org/sequence/id/$i?content-type=text/plain`; done > 106Edger_Transcritos.fa

