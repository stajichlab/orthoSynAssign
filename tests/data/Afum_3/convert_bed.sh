 zcat Afum_3/gff3/FungiDB-68_AfumigatusAf293.gff3.gz | grep -v ^# | grep -P protein_coding_gene | cut -f1,4,5,9 | perl -p -e 's/ID=([^;]+);.+/$1/' > Afum_3/bed/FungiDB-68_AfumigatusAf293.bed
zcat FungiDB-68_AnovofumigatusIBT16806.gff3.gz |  grep -v ^# | grep -P protein_coding_gene | cut -f1,4,5,9 | perl -p -e 's/ID=([^;]+);.+/$1/' > ../bed/FungiDB-68_AnovofumigatusIBT16806.bed
