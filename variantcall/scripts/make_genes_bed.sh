pushd genome
grep -P "\tgene\t" FungiDB-39_AfumigatusAf293.gff | cut -f 1,4,5,9 | perl -p -e 's/ID=([^;]+);.+/$1/' | sort -k1,1 -k2,2n > Af293.genes.bed
popd
